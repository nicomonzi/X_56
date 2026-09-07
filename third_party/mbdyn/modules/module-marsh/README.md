# module-marsh

A module enabling MBDyn as a **Flight Dynamics Module** node for
[MARSH](https://marsh-sim.github.io/)
(Modular Aircraft Research Simulator Handler).

MARSH is a modular flight-simulator framework whose nodes (flight model,
pilot controls, visualisation, instruments, motion platform, ...) exchange
[MAVLink](https://mavlink.io/) messages over UDP, routed by the `MARSH
Manager`. This module lets an MBDyn model take the flight-model role: it
publishes the aircraft state computed by MBDyn and acquires pilot inputs,
making them available to the rest of the model as ordinary drive-caller
values.

This README is aimed at developers of the module. For end-user input syntax
see the MBDyn input manual, section *Modules / Element Modules /
Module-marsh* (`manual/input/module-marsh.tex`).

## Quick start

```
./configure --with-module="marsh" --enable-runtime-loading ...
make && make install
```

Minimal deck (see `marsh_test.mbd` for a complete, runnable one):

```
module load: "libmodule-marsh";

air properties: std, SI, null;

aircraft instruments: 1, 1,
    orientation, flight mechanics,
    initial latitude, 45.5*deg2rad,
    initial longitude, 9.15*deg2rad;

user defined: 3, marsh,
    manager address, "127.0.0.1",
    manager port, 24400,
    send, "SIM_STATE",
        aircraft instruments, 1,
    send, "MOTION_CUE_EXTRA",
        aircraft instruments, 1,
    subscribe, "MANUAL_CONTROL";

# use a received pilot input anywhere in the model
drive caller: 100,
    element, 3, loadable, string, "MANUAL_CONTROL.x", direct;
```

## Files

| File | Contents |
| --- | --- |
| `module-marsh.{h,cc}` | `ModuleMarsh`, the `UserDefinedElem`: input parsing, UDP socket, per-step send/receive loop, heartbeat scheduling, private-data fan-out, `module_init` |
| `mavlink_pubsub.h` | The `MavlinkProducer` / `MavlinkConsumer` interfaces and the two registry lookup functions |
| `mavlink_producers.{h,cc}` | Outbound messages (`SIM_STATE`, `MOTION_CUE_EXTRA`), the `MarshMotionSource` state sources, the per-field override machinery, and the producer registry |
| `mavlink_consumers.{h,cc}` | Inbound messages (`MANUAL_CONTROL`) plus their registry table |
| `mavlink/` | Vendored, generated MAVLink v2 C headers (see below) |
| `Makefile.inc` | Extra objects and compiler flags for the generic module build |
| `marsh_test.mbd` | Regression test used by the CI `mbdyn-modules-test-job` |

## Architecture

### MAVLink Producers / Consumers

`ModuleMarsh` owns a UDP socket and two vectors of handlers:

- **producers** — `MavlinkProducer::Pack()` fills one outbound
  `mavlink_message_t` from live MBDyn state;
- **consumers** — `MavlinkConsumer::Decode()` unpacks one inbound message
  into member fields, and optionally exposes them as private data.

Which handlers exist in a given element instance is decided at parse time
by the `send` / `subscribe` clauses in the input file, each naming a MAVLink
message. The name is looked up in a `std::map<std::string, Factory>` table
(`ProducerRegistry` / `ConsumerRegistry`), and the selected factory picks
any message-specific arguments from the parser. `HEARTBEAT` is added as
a protocol-level feature, always sent and always processed.

This keeps the element itself message-agnostic — it never mentions
`SIM_STATE` or `MANUAL_CONTROL` — so adding a message touches only one
`.cc` file.

### Message sources

Where an outbound message takes its state from is a `MarshMotionSource`,
which reports everything in MAVLink conventions (body axes
x-forward/y-right/z-down, SI units, NED for world components) and declares
what it cannot provide:

- **`InstrumentsMotionSource`** wraps an `aircraft instruments` element and
  reads its private data. This is the source that fills in the whole of
  `SIM_STATE`: the elements already resolves the aircraft frame including
  its offset with respect to the node, the flight-mechanics angles, and the
  geodetic position. Private-data indices are resolved once, at parse time,
  via `iGetPrivDataIdx()`.
- **`NodeMotionSource`** binds a `StructNode` plus an optional offset and
  rotation, in `module-imu` style. For the cases that need no
  instruments element. It has no geodetic position and no flight-mechanics
  Euler angles, so those fields stay null unless set explicitly.

A source is optional: a message can be assembled out of drive callers alone
(see the manual entry and the next section).

### Per-field overrides

Any float field of a message can be driven by an arbitrary `ScalarValue`
(`drive, <DriveCaller>` or `node dof, <ScalarDof>`) with a
`field, "<name>", <value>` clause, reusing the core `ReadScalarValue()`
helper from `mbdyn/base/scalarvalue.h` — the same one module-flightgear
uses. Overrides are applied after the source has filled the message, so
they win.

The mechanism is a table of pointers to member, one per message:

```c++
static const MavlinkFieldDesc<mavlink_sim_state_t> SimStateFields[] = {
        { "q1", &mavlink_sim_state_t::q1 },
        ...
```

which is why producers build a `mavlink_<msg>_t` struct and call
`mavlink_msg_<msg>_encode()`, rather than the positional `_pack()` form.

This exists because the interface is useful beyond aircraft dynamics: a
field can carry a buffet, a vibration, a simulated sensor fault, or any
other effect the flight model does not itself describe.

### Per-step cycle

Everything happens in `ModuleMarsh::AfterConvergence()`, once per converged
time step:

1. `ReceivePending()` — drain the socket with non-blocking `recvfrom()`
   until `EAGAIN`, feeding every byte through `mavlink_parse_char()`. This
   handles several messages coalesced into one datagram, and only completed,
   CRC-checked messages are dispatched. Each decoded message goes to
   `Dispatch()`, which either updates the manager-seen flag (`HEARTBEAT`) or
   is offered to every consumer whose `GetMsgId()` matches.
2. `SendHeartbeat()` — sends a `HEARTBEAT` if at least
   `1 / heartbeat rate` of *simulation* time has elapsed since the last one.
3. `SendProducers()` — packs and sends every registered producer.

There is no timer thread and no OS timer. MBDyn's real-time POSIX mode
(`RTPOSIXSolver`) paces the *solver's* step loop with `clock_nanosleep()`,
and elements simply see `AfterConvergence()` once per step.
Comparing `m_pDM->dGetTime()` against a stored timestamp therefore
behaves correctly both under real-time pacing and in batch runs.

### Coordinate frames, units and the gravity convention

MBDyn's flight-mechanics "world" is North-West-Up, while MAVLink's is
North-East-Down, hence the sign flips in `GetVelNED()`. Body axes agree
between the two (x forward, y right, z down), so body-frame quantities pass
straight through. Geodetic coordinates are radians in MBDyn and degrees in
MAVLink, converted when packing, along with the `lat_int`/`lon_int`
high-precision fields in degE7.

MAVLink's `SIM_STATE` and MARSH's `MOTION_CUE_EXTRA` report **specific force**
 — what an accelerometer reads, i.e. gravity included — not acceleration.
A body at rest reports a vector of magnitude *g*.
`GetSpecificForce()` therefore computes

```
Rᵀ · (XPP − gravity)
```

with `gravity` obtained from `GravityOwner::bGetGravity(X, Acc)`, which
leaves its output untouched (so, zero) when the deck defines no `gravity`
element. The reference MARSH node `mps-adapter` performs the mirror-image
subtraction when *consuming* these messages, so the two agree.

Note that specific force is computed here rather than being read from the
`aircraft instruments` element, because gravity is not reachable from that
element: `DataManager` wires the gravity pointer only into element types
whose `bGeneratesInertiaForces()` is set — `BODY`, `SOLID`, `JOINT`,
`LOADABLE` and `PLATE` (`mbdyn/base/elman.cc`) — and `aircraft instruments`
is `AERODYNAMIC`. This module, being a `UserDefinedElem`, is `LOADABLE`, so it
does get gravity.

### Private data

`ModuleMarsh` does not own any private-data names itself; it concatenates
the ranges published by its consumers. `iGetPrivDataIdx()` offers the
requested string to each consumer in turn and biases the returned local
index by the running offset; `dGetPrivData()` walks the same offsets in
reverse via `GetPrivDataOwner()`.

Consumers are passed the **full dotted string** (`"MANUAL_CONTROL.x"`).
Keeping the message name in the key is what makes the flat namespace
unambiguous once several consumers are registered.

## Changes made to `aircraft instruments`

Supporting this module required extending `aircraft instruments`
(`mbdyn/aero/instruments.{h,cc}`), which is the natural source of flight
state and already provided most of what `SIM_STATE` needs:

- **Fixed** `body_axb` / `body_ayb` / `body_azb` and `body_pd` / `body_qd` /
  `body_rd`: they are documented as being in the body frame, but stored the
  components of `GetXPPCurr()` / `GetWPCurr()` as they came, i.e. in the
  global frame. They are now rotated into the aircraft frame.
- **Added** `e0` / `e1` / `e2` / `e3`, the attitude of the aircraft frame as
  Euler parameters. Unlike `attitude` / `bank` / `heading`, which are built
  from `asin()` of single matrix entries, these are unambiguous — and they
  expose the aircraft frame *including* the element's own orientation
  offset, which is otherwise not reachable from outside.
- **Added** `vnorth` / `veast` / `vdown`, proper NED components of the node
  velocity. Note these are ground velocity, whereas `verticalspeed`
  subtracts the airstream.

## Adding a MAVLink message

To publish a new outbound message:

1. Derive from `MavlinkProducer` in `mavlink_producers.h`, holding a
   `MarshMotionSource *` and a vector of overrides if it should support them.
2. Add a `MavlinkFieldDesc<mavlink_<msg>_t>` table naming its float fields.
3. Implement the constructor — typically just `ReadMarshMotionSource()` plus
   `ReadMavlinkFieldOverrides()` — and `Pack()`, filling a
   `mavlink_<msg>_t`, applying the overrides, then calling
   `mavlink_msg_<msg>_encode()`.
4. Add one line to `ProducerRegistry` in `mavlink_producers.cc`.

For a new inbound message the shape is the same against `MavlinkConsumer`,
`ConsumerRegistry` in `mavlink_consumers.cc`, plus `GetMsgId()`,
`Decode()` (via `mavlink_msg_<name>_decode()`) and, if the values should
reach the model, the three private-data methods.

## Vendored MAVLink headers

`mavlink/` holds generated MAVLink v2 C headers. Four dialect directories
are required, because they form an include chain:

```
marsh/marsh.h → common/common.h → standard/standard.h → minimal/minimal.h
```

`HEARTBEAT` lives in `minimal`. The MARSH-specific messages
(`MOTION_CUE_EXTRA`, `MOTION_PLATFORM_STATE`, `CONTROL_LOADING_AXIS`, ...)
and the `MARSH_COMP_ID_*` component-id enum are in `marsh/`.

Include the chain through its dialect entry point, *after* `dataman.h`:

```c
#include "mavlink/marsh/mavlink.h"
```

To refresh against upstream, regenerate with `pymavlink`'s `mavgen` from
[`marsh-sim/mavlink`](https://github.com/marsh-sim/mavlink) (branch
`dialect`, `message_definitions/v1.0/all.xml`, wire protocol 2.0); the
`scripts/update_mavlink.py` helper in `mps-adapter` does exactly this and
can be pointed here. The dialect XML is the authoritative field-by-field
reference, mirrored as generated docs at
<https://marsh-sim.github.io/mavlink/>.

The generated code uses packed structs, so `Makefile.inc` adds
`-Wno-address-of-packed-member`.

## MARSH conventions

- **Component ID** identifies the node's *role*; this module defaults to
  `MARSH_COMP_ID_FLIGHT_MODEL` (26). The Manager forwards to a component
  only the messages relevant to its ID.
- **System ID** identifies the simulation *session*: all nodes taking part
  in one simulation should share it (default 1).
- Only the first component announcing a given (system id, component id)
  pair has its messages forwarded; later ones are shadowed by the Manager.
- The Manager listens on UDP port **24400** by default.
- `HEARTBEAT` is expected at 1 Hz.

## Testing

`marsh_test.mbd` is a deliberately trivial deck — a clamped rigid body under
gravity — chosen so the published values are constants that are easy to
eyeball: the body never moves, so the specific force is steady and the body
rates are zero. It covers both message sources (`SIM_STATE` from the
instruments element, `MOTION_CUE_EXTRA` from the node) and one drive-caller
field override. A `genel clamp` on an abstract node pulls
`MANUAL_CONTROL.x` into the output files.

At the bottom of `marsh_test.mbd`, the `## @MBDYN_SIMPLE_TESTSUITE_EXIT_STATUS@`
trailer that `testsuite/simple_testsuite.sh` uses as its pass/fail baseline; see
the related [CI wiki page](https://public.gitlab.polimi.it/DAER/mbdyn/-/wikis/CI).

For a round-trip check against a real peer, run the deck alongside the MARSH
Manager, or against a `pymavlink` script that listens for the
outbound messages and injects a `MANUAL_CONTROL`.

## Current limitations

- **No manager watchdog.** A `HEARTBEAT` from the Manager sets a flag and a
  timestamp, and the flag is reported in the output file, but nothing acts
  on the Manager going silent. Should treat 5 s without a heartbeat
  as a disconnection.
- **Producers send every step.** There is no per-message rate limit, so a
  small time step means a high datagram rate. A `send every` style option,
  mirroring the `output every` of MBDyn's own stream output elements, would
  be the natural fix.
- **`SIM_STATE` standard deviations are null**, as MBDyn models no sensor
  noise; they can be set with `field` clauses.
- **`MANUAL_CONTROL` axes are raw.** Values reach the model as integers
  in [-1000, 1000]; scaling is left to the drive caller wrapping the
  private data.
- Only float message fields can be overridden with `field` clauses; integer
  fields such as `time_boot_ms` cannot.
