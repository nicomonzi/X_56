# X-56A MBDyn-DUST final BFF diagnostic

This directory is the Git-portable CHECK/SMOKE/PRODUCTION architecture for one
diagnostic at 60.8421 m/s. Production is deliberately blocked by the runner
until the scientific blockers listed below are resolved. Never use SMOKE output
for aerodynamic, trim, damping, or flutter conclusions.

## New clone / university server

```bash
git clone <repository-url>
cd TESI/BFF_DUST_55
cp config/machine.env.example config/machine.env
# Edit only executable locations and MBDYN_PYTHON_PATH.
python3 run_case.py --check-only
python3 run_case.py --smoke --threads 8
```

After all production blockers have been resolved and the audit report has been
regenerated:

```bash
python3 run_case.py --production --threads 12
```

Do **not** run production on the development machine unless explicitly
requested. No package installation, compilation, MPI launch, sweep, or time-step
sensitivity is performed by the launcher.

## Execution levels

### CHECK

`python3 run_case.py --check-only`

CHECK verifies executables, preCICE support/API, Python requirements, required
files, the COARSE mesh, coupling-node count, portable paths, writable output,
and lets MBDyn and DUST parse to their coupling waits. It advances zero physical
time steps.

### SMOKE

`python3 run_case.py --smoke --threads 8`

SMOKE uses 60.8421 m/s, the 900-panel reconstructed COARSE mesh, `dt=0.002 s`,
and `t_end=0.10 s`: exactly 50 physical windows. It has no complete trim, RAP,
open-loop identification, Matrix Pencil analysis, controller tuning, or sweep.
Every result is labelled `SMOKE TEST - NOT PHYSICALLY VALID`. The authorized
local smoke completed all 50 converged windows, generated a continuous wake,
and terminated both solvers normally in about 127 s wall time using 8 threads.
Server timing may differ.

### PRODUCTION

`python3 run_case.py --production --threads 12`

The prepared candidate uses the reconstructed 4284-panel FINE mesh,
`dt=0.002 s`, `t_end=9.50 s`, and 4750 physical windows. The state sequence is
CAPTURE_TRIM → READY → RAP → OPEN_LOOP → RECOVERY → END. DUST and its wake are
never restarted between states. The STRIP/C81 trim is only an initial guess;
the active longitudinal SAS acquires/holds the coupled DUST state before the
RAP. During OPEN_LOOP the sample-and-hold disconnects feedback and freezes all
ten actual surface commands. Feedback reconnects only in RECOVERY.

The estimated storage is approximately 5–40 GiB, mainly dependent on actual
vortex-particle population and DUST HDF5 compression. There are 474 requested
DUST output instants at 0.02 s (`output_start=F`); MBDyn/coupling scalar histories remain at
0.002 s. Keep at least 50 GiB free until a short server run measures the actual
rate.

## CPU and NUMA policy

`--threads` accepts only 8, 12, or 16 and defaults to 12. Before launch:

```text
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=<requested>
OMP_PLACES=cores
OMP_PROC_BIND=close
```

The runner prints hostname, date, uptime, available RAM, load average, and all
thread settings. It never increases the requested count and does not use MPI.

## Modal basis conflict

The available FEM has six near-zero rigid modes followed by elastic FEM modes:

- FEM 7: 3.21713 Hz, identified as the first symmetric wing-bending mode;
- FEM 18: 18.55899 Hz;
- FEM 24: 25.02676 Hz.

Consequently, FEM 7–18 contains **12 elastic modes**, while 18 elastic modes
would require FEM 7–24 and would violate the requested 18.554 Hz upper limit.
The provisional model respects the frequency cap and retains FEM 7–18.
Production remains blocked until the intended interpretation is confirmed.
The complete provisional list is in `reference/modal_basis.csv` and is printed
by CHECK.

## Mesh audit and current blockers

`reports/fine_mesh_audit.json` and `reports/fine_mesh_audit.vtu` are generated
with:

```bash
python3 tools/audit_mesh.py
```

The FINE mesh is labelled **production diagnostic mesh with residual
convergence uncertainty**, never mesh-converged. Current blockers are:

1. the available parametric topology has no winglet component;
2. the original mesh-study input files were absent, so 900/2280/4284 levels
   were reconstructed without changing the available planform;
3. not all hinge endpoints coincide with spanwise mesh lines;
4. the supplied MEDIUM→FINE differences remain 3.16% in Fz, 14.59% in My, and
   26.33% in spanwise loading;
5. the modal count/frequency requirement is internally inconsistent;
6. the rigid DUST control-effectiveness diagnostic has not been executed.

The rigid diagnostic design is in `reference/control_effectiveness_plan.csv`.
It requires baseline, theta ±1°, and symmetric body-flap ±1° COARSE cases. Put
their integrated values in `results/control_effectiveness/loads.csv` and run:

```bash
python3 tools/analyze_control_effectiveness.py
```

Do not compensate a failed sign/axis/reference check with stronger SAS gains.

## Wake settings

| Setting | SMOKE | PRODUCTION candidate |
|---|---:|---:|
| `n_wake_panels` | 4 | 40 |
| particle box [in] | 1500 × 1000 × 1000 | 26000 × 3000 × 3000 |
| particle limit | 50,000 | 1,500,000 |
| update frequency | every step | every step |
| FMM | enabled | enabled |
| current particle count | unavailable before run | unavailable before run |

The production downstream extent contains approximately 9.5 s of convection
at the specified velocity. Particle count must be read from the first server
result; the configured value is a limit, not a predicted population.

## preCICE portability

Both the tested legacy 2.x configuration and a ported 3.x candidate are kept.
The runner selects by `PRECICE_MAJOR` or Python API detection. The university
server must validate the 3.2 XML and DUST build during CHECK before production.
The 3.x adapter uses `precice.Participant`, named mesh/data access, and the v3
checkpoint API. No v3.2 coupled run has been performed on this development
machine.

## Output required for identification

`diagnostics.csv` consolidates each accepted MBDyn text state with the preCICE
diagnostics at the 0.002 s rate: theta, q, Vz, Z, SWB1 coordinate/velocity,
symmetric/antisymmetric tips, all ten surface angles, global Fz/My, state, and
iteration count. The active development build leaves the NetCDF time dimension
empty during iterative coupling, so `tools/extract_mbdyn_text.py` deliberately
uses `.mov/.mod/.jnt` as the authoritative structural source. Future runs form
global aerodynamic moments about the configured CG; the already completed
smoke is explicitly marked as containing the legacy sum of nodal moments only.
DUST does not expose a live particle count through this coupling API, so that
column remains blank and must be populated from DUST postprocessing on the
server. Production postprocessing must verify surface continuity and zero
feedback during OPEN_LOOP before Matrix Pencil results are accepted.

ParaView time collections are written under `output/<level>/paraview/`. Open
`x56_coupled.pvd` for the deformed aerodynamic surface,
`x56_coupled_wpan.pvd` for wake panels, and `x56_coupled_wpart.pvd` for vortex
particles. The static MBDyn mapping is `reports/fine_coupling_nodes.vtu`.
