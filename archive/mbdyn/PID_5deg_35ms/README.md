# X-56 impulsive PID tests

This folder is a self-contained copy of `MBDYN/PID` with two isolated-axis tests:

- `main_roll.mbd`: only roll is free; left and right surfaces move differentially.
- `main_pitch.mbd`: only pitch is free; all movable surfaces move symmetrically.

The roll case uses a 10° command, 40 m/s airspeed, and ±30° surface limits.
The pitch case uses a 5° command, 30 m/s airspeed, and ±10° surface limits.
Both targets are applied as instantaneous steps at 2.0 s, without angular-rate
constraints. Both cases use a 0.02 s time step and a 10 s final time.
The pin joint acts directly on the center-of-gravity node. Gravity is disabled
only in the pitch-isolation case, where it would otherwise create a permanent
moment about the single free axis. The pitch controller combines an angle PI
(`Kp=0.32`, `Ki=0.68`) with a separate rigid-body pitch-rate damper
(`Kp=0.2`). Their limits are ±7° and ±3°, so the combined command remains
inside the ±10° actuator envelope.
The underdamped roll tuning uses `Kp=3.6`, `Ki=3.1`, and `Kd=1.3`, also
with `Kn=15`. It produces a moderate recovered overshoot; the differential
surfaces touch the ±30° limits only during the initial command transient.

Generated dashboards:

- `dashboard_roll_pid_10deg_40ms.png`
- `dashboard_pitch_pid_5deg_30ms.png`

The corresponding simulation results are `output/roll_10deg_40ms.nc` and
`output/pitch_5deg_30ms.nc`.

Run both simulations and generate the English dashboards with:

```bash
python3 run_pid_tests.py
```

To regenerate plots from existing NetCDF files:

```bash
python3 plot_pid_results.py --case all
```
