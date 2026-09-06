#!/usr/bin/env python3
"""Run one bounded NASA-style X-56 full-surface-hold SAS-off observation."""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / "main_bff_open_loop.mbd"
DEFAULT_OUTPUT = Path("/mnt/c/Users/Utente/Desktop/BFF_open_loop")
OUTPUT_ENV = "BFF_OPEN_LOOP_OUTPUT_DIR"
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")
CASE_SUFFIXES = (".mbd", ".nc", ".out", ".log", ".stdout", ".mov", ".usr", ".aer")

# NASTRAN-calibrated generalized correction for FEM mode 7.  K7 and C7 are
# not flight-control gains: they compensate the missing 3-D DLM/RFA terms in
# the sectional Theodorsen model.  Values between knots are interpolated
# linearly; the sweep is intentionally restricted to this validated range.
DLM_ROM_SCHEDULE = (
    # TAS [m/s], K7 [1/s^2], C7 [1/s], refined point-7 frequency [Hz]
    (50.0, 40.0000, 0.0000, 2.7547653),
    (52.5, 50.0000, 0.0000, 2.6673810),
    (55.0, 72.0000, -0.2000, 2.5520311),
    (57.5, 109.0000, -0.8800, 2.3787639),
    (60.0, 142.0000, 2.5000, 2.1063462),
    (62.5, 155.0000, 6.0000, 2.0605652),
    (65.0, 152.0000, 8.4700, 2.0551998),
    (67.5, 145.0000, 11.2000, 2.0584934),
    (70.0, 146.0000, 12.8000, 2.0643811),
)

# Closed-loop flutter-suppression authority used only to enter/recover the
# point.  It is multiplied by the SAS gate and is exactly zero in open loop.
# Above 65 m/s the released plant needs progressively more mode-7 damping.
SAS_MODAL_DAMPING_SCHEDULE = (
    (50.0, 16.0),
    (65.0, 16.0),
    (67.5, 32.0),
    (70.0, 50.0),
)


def _interpolate_schedule(velocity: float, column: int) -> float:
    if not DLM_ROM_SCHEDULE[0][0] <= velocity <= DLM_ROM_SCHEDULE[-1][0]:
        raise ValueError("la correzione SOL 145 e' validata solo tra 50 e 70 m/s")
    for lower, upper in zip(DLM_ROM_SCHEDULE, DLM_ROM_SCHEDULE[1:]):
        if lower[0] <= velocity <= upper[0]:
            fraction = (velocity - lower[0]) / (upper[0] - lower[0])
            return lower[column] + fraction * (upper[column] - lower[column])
    return DLM_ROM_SCHEDULE[-1][column]


def dlm_rom_values(velocity: float) -> tuple[float, float, float, float]:
    k7 = _interpolate_schedule(velocity, 1)
    c7 = _interpolate_schedule(velocity, 2)
    # Quadratic fit of the settled mode-7 coordinate from the current
    # Theodorsen/3/4-chord model; it only removes the static trim component.
    q7_equilibrium = (
        -7.10495281e-05 * velocity * velocity
        + 1.10266063e-02 * velocity
        - 8.38346238e-01
    )
    rap_frequency = _interpolate_schedule(velocity, 3)
    return k7, c7, q7_equilibrium, rap_frequency


def sas_modal_damping(velocity: float) -> float:
    if not SAS_MODAL_DAMPING_SCHEDULE[0][0] <= velocity <= SAS_MODAL_DAMPING_SCHEDULE[-1][0]:
        raise ValueError("lo scheduling SAS modale e' definito solo tra 50 e 70 m/s")
    for lower, upper in zip(SAS_MODAL_DAMPING_SCHEDULE, SAS_MODAL_DAMPING_SCHEDULE[1:]):
        if lower[0] <= velocity <= upper[0]:
            fraction = (velocity - lower[0]) / (upper[0] - lower[0])
            return lower[1] + fraction * (upper[1] - lower[1])
    return SAS_MODAL_DAMPING_SCHEDULE[-1][1]


def validate_release_state(nc_path: Path, rendered: str) -> None:
    """Reject a BFF run that did not enter the test point quasi-steadily."""
    constants = {
        name: float(value)
        for name, value in re.findall(
            r"(?m)^set:\s*const\s+real\s+(SAS_OFF_START|TIME_STEP)\s*=\s*([-+0-9.eE]+)\s*;",
            rendered,
        )
    }
    with Dataset(nc_path) as data:
        time = np.asarray(data["time"][:]).squeeze()
        velocity = np.asarray(data["node.struct.990000.XP"][:])
        omega = np.asarray(data["node.struct.990000.Omega"][:])
    release = int(np.searchsorted(time, constants["SAS_OFF_START"]))
    if release >= len(time):
        raise RuntimeError(f"{nc_path.name}: run interrotta prima del rilascio SAS")
    window = (time >= constants["SAS_OFF_START"] - 1.0) & (time < constants["SAS_OFF_START"])
    vz = velocity[:, 2] * 0.0254
    q = omega[:, 1] * 180.0 / math.pi
    p = omega[:, 0] * 180.0 / math.pi
    checks = {
        "abs_Vz_release_mps": abs(float(vz[release])),
        "abs_q_release_deg_s": abs(float(q[release])),
        "abs_p_release_deg_s": abs(float(p[release])),
        "std_Vz_last_1s_mps": float(np.std(vz[window])),
        "std_q_last_1s_deg_s": float(np.std(q[window])),
    }
    limits = {
        "abs_Vz_release_mps": 0.05,
        "abs_q_release_deg_s": 0.10,
        "abs_p_release_deg_s": 0.20,
        "std_Vz_last_1s_mps": 0.05,
        "std_q_last_1s_deg_s": 0.20,
    }
    failures = [f"{name}={checks[name]:.6g}>{limit:.6g}" for name, limit in limits.items() if checks[name] > limit]
    if failures:
        raise RuntimeError(f"{nc_path.name}: rilascio non trimmato: " + ", ".join(failures))
    print(
        "[trim] release valido: "
        f"Vz={vz[release]:+.5f} m/s, q={q[release]:+.5f} deg/s, p={p[release]:+.5f} deg/s"
    )


def case_stem(velocity: float, settling_probe: float | None = None) -> str:
    stem = f"NASA_OL_V_{velocity:08.4f}".replace(".", "p")
    if settling_probe is not None:
        stem += f"_SETTLING_{settling_probe:05.2f}s".replace(".", "p")
    return stem


def render_case(velocity: float, settling_probe: float | None = None) -> str:
    if not math.isfinite(velocity) or velocity <= 0.0:
        raise ValueError("V_INF deve essere positiva e finita")
    text, count = re.subn(
        r"(?m)^set:\s*const\s+real\s+V_INF\s*=\s*[^;]+;",
        f"set: const real V_INF = {velocity:.8f};",
        MODEL.read_text(),
    )
    if count != 1:
        raise RuntimeError(f"attesa una definizione di V_INF, trovate {count}")
    k7, c7, q7_equilibrium, rap_frequency = dlm_rom_values(velocity)
    replacements = {
        "NASTRAN_DLM_ROM_K7": k7,
        "NASTRAN_DLM_ROM_C7": c7,
        "NASTRAN_DLM_ROM_Q7_EQ": q7_equilibrium,
        "SAS_MODAL_DAMPING": sas_modal_damping(velocity),
    }
    setconst = (ROOT / "INCLUDE/setconst.mbd").read_text()
    for name, value in replacements.items():
        setconst, replaced = re.subn(
            rf"(?m)^set:\s*const\s+real\s+{name}\s*=\s*[^;]+;",
            f"set: const real {name} = {value:.10e};",
            setconst,
        )
        if replaced != 1:
            raise RuntimeError(f"attesa una definizione di {name}, trovate {replaced}")
    text, replaced = re.subn(
        r"(?m)^set:\s*const\s+real\s+BFF_RAP_TARGET_FREQUENCY\s*=\s*[^;]+;",
        f"set: const real BFF_RAP_TARGET_FREQUENCY = {rap_frequency:.10e};",
        text,
    )
    if replaced != 1:
        raise RuntimeError("definizione BFF_RAP_TARGET_FREQUENCY non univoca")
    if settling_probe is not None:
        if not math.isfinite(settling_probe) or settling_probe <= 0.0:
            raise ValueError("la durata del settling probe deve essere positiva e finita")
        # Keep the SAS continuously engaged.  The open-loop event is moved
        # beyond the final time, so this diagnostic cannot be mistaken for a
        # BFF identification run.
        text, replaced_final = re.subn(
            r"(?m)^set:\s*const\s+real\s+FINAL_TIME\s*=\s*[^;]+;",
            f"set: const real FINAL_TIME = {settling_probe:.8f};",
            text,
        )
        text, replaced_off = re.subn(
            r"(?m)^set:\s*const\s+real\s+SAS_OFF_START\s*=\s*[^;]+;",
            f"set: const real SAS_OFF_START = {settling_probe + 1.0:.8f};",
            text,
        )
        if replaced_final != 1 or replaced_off != 1:
            raise RuntimeError("costanti temporali del settling probe non univoche")
    required = (
        "set: const real TIME_STEP = 0.01;",
        "sample and hold",
        '"((Time<SAS_OFF_START)||(Time>=SAS_ON_START))"',
        "SAS_OFF_DURATION = 2.05",
        "BFF_RAP_DURATION = 0.742/BFF_RAP_TARGET_FREQUENCY",
        "force: NASTRAN_DLM_ROM_FORCE, modal, MODAL_JOINT",
        "force: SAS_MODAL_DAMPER_FORCE, modal, MODAL_JOINT",
        "LIFT_HOLD_GAIN =",
        "reference, LAT_LIMIT_DRIVE",
        "reference, BFL_HOLD_DRIVE",
        "reference, BFR_HOLD_DRIVE",
        "reference, WF1L_HOLD_DRIVE",
        "reference, WF1R_HOLD_DRIVE",
        "reference, WF2L_HOLD_DRIVE",
        "reference, WF2R_HOLD_DRIVE",
        "reference, WF3L_HOLD_DRIVE",
        "reference, WF3R_HOLD_DRIVE",
        "reference, WF4L_COMMAND_DRIVE",
        "reference, WF4R_COMMAND_DRIVE",
    )
    combined = (
        text
        + setconst
        + (ROOT / "INCLUDE/control_surfaces.mbd").read_text()
    )
    if not all(item in combined for item in required):
        raise RuntimeError("invarianti open-loop/dt non soddisfatte")
    forbidden = ("NOTCH_", "BFF_BS", "SAFETY_ACT_DRIVE", "maneuver", "MANEUVER")
    if any(item in combined for item in forbidden):
        raise RuntimeError("il caso contiene notch, safety o logica di manovra non ammessi")
    aerobody = (ROOT / "INCLUDE/aerobody.mbd").read_text()
    if aerobody.count("theodorsen, c81") != 58:
        raise RuntimeError("tutti i 58 elementi devono usare Theodorsen/C81")
    if aerobody.count("-0.5*(") != 116:
        raise RuntimeError("collocazione a 3/4 di corda incompleta nei 58 elementi")
    rendered = text.replace('"./INCLUDE/', f'"{ROOT / "INCLUDE"}/')
    # The rendered case must be self-contained with respect to scheduled ROM
    # values even though the remaining include files stay shared.
    rendered_setconst = ROOT / "INCLUDE/setconst.mbd"
    rendered = rendered.replace(
        f'include: "{rendered_setconst}";',
        setconst,
    )
    return rendered


def run_case(
    velocity: float,
    output: Path,
    overwrite: bool = False,
    settling_probe: float | None = None,
) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    stem = case_stem(velocity, settling_probe)
    prefix = output / stem
    nc_path = prefix.with_suffix(".nc")
    rendered = render_case(velocity, settling_probe)
    input_path = prefix.with_suffix(".mbd")
    if nc_path.exists() and not overwrite:
        if input_path.exists() and input_path.read_text() == rendered:
            print(f"[reuse] {stem}")
            return nc_path
        raise RuntimeError(f"{stem}: output esistente non coerente; usare --overwrite")
    if overwrite:
        for suffix in CASE_SUFFIXES:
            path = output / f"{stem}{suffix}"
            if path.is_file():
                path.unlink()
    input_path.write_text(rendered)
    command = [str(MBDYN), "-s", "-f", str(input_path), "-o", str(prefix)]
    protocol = "closed-loop settling probe" if settling_probe is not None else "bounded SAS-off"
    print(f"[run] {stem}: V={velocity:.4f} m/s, {protocol}")
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    prefix.with_suffix(".stdout").write_text(result.stdout + result.stderr)
    if result.returncode or not nc_path.exists():
        raise RuntimeError(f"MBDyn fallito ({result.returncode}); vedere {prefix.with_suffix('.stdout')}")
    if settling_probe is None:
        validate_release_state(nc_path, rendered)
    print(f"[done] {nc_path}")
    return nc_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("velocity", nargs="?", type=float, default=60.8421, help="V_INF [m/s]")
    parser.add_argument("--output", type=Path, default=Path(os.environ.get(OUTPUT_ENV, DEFAULT_OUTPUT)))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--settling-probe", type=float, metavar="SECONDS",
        help="run only a closed-loop settling diagnostic; no SAS-off event",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_case(args.velocity, args.output, args.overwrite, args.settling_probe)


if __name__ == "__main__":
    main()
