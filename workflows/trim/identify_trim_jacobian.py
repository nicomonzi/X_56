#!/usr/bin/env python3
"""Identifica il Jacobiano 2x2 del trim con perturbazioni finite centrate.

La configurazione aerodinamica e' quella di main_trim.mbd: soltanto BFL/BFR
ricevono la deflessione simmetrica, mentre tutti i wing flap restano a zero.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from run_trim_sweep import (
    ROOT,
    completed_case,
    locate_mbdyn,
    replace_real_constant,
)


DEFAULT_OUTPUT = ROOT / "jacobian_calibration_bfl_bfr"
REFERENCE_SPEED = 63.0
REFERENCE_DENSITY = 1.146e-7
SETTLED_START = 22.0
RAD_TO_DEG = 180.0 / math.pi


@dataclass(frozen=True)
class Perturbation:
    name: str
    pitch_deg: float
    elevator_deg: float


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "main_trim.mbd")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mbdyn", type=Path, default=None)
    parser.add_argument("--pitch-step", type=float, default=0.25, help="[deg]")
    parser.add_argument("--elevator-step", type=float, default=0.50, help="[deg]")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def calibration_source(
    template: str, pitch_deg: float, elevator_deg: float
) -> str:
    source = replace_real_constant(template, "V_INF", REFERENCE_SPEED)
    source = replace_real_constant(source, "RHO_AIR", REFERENCE_DENSITY)
    source = replace_real_constant(source, "TRIM_RESPONSE_RATE", 0.0)
    source = replace_real_constant(
        source, "TRIM_PITCH_INITIAL", math.radians(pitch_deg)
    )
    return replace_real_constant(
        source, "TRIM_ELEVATOR_INITIAL", math.radians(elevator_deg)
    )


def launch_case(
    mbdyn: Path,
    template: str,
    output_directory: Path,
    perturbation: Perturbation,
    force: bool,
) -> tuple[str, int, float]:
    case_directory = output_directory / perturbation.name
    case_directory.mkdir(parents=True, exist_ok=True)
    if not force and completed_case(case_directory):
        return perturbation.name, 0, 0.0

    archive_input = case_directory / "case_input.mbd"
    archive_input.write_text(
        calibration_source(
            template, perturbation.pitch_deg, perturbation.elevator_deg
        ),
        encoding="utf-8",
    )
    temporary_input = ROOT / f".__jacobian_{perturbation.name}.mbd"
    temporary_input.write_text(archive_input.read_text(encoding="utf-8"), encoding="utf-8")
    command = [str(mbdyn), "-f", temporary_input.name, "-o", str(case_directory / "result")]
    started = time.monotonic()
    try:
        process = subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    finally:
        temporary_input.unlink(missing_ok=True)
    elapsed = time.monotonic() - started
    (case_directory / "console.log").write_text(
        process.stdout, encoding="utf-8", errors="replace"
    )
    return perturbation.name, process.returncode, elapsed


def variable(nc: Dataset, name: str) -> np.ndarray:
    return np.ma.filled(nc.variables[name][:], np.nan).astype(float)


def read_result(
    output_directory: Path, perturbation: Perturbation
) -> dict[str, float | str]:
    with Dataset(output_directory / perturbation.name / "result.nc") as nc:
        time_values = variable(nc, "time")
        pitch = variable(nc, "node.struct.990000.Phi")[:, 1]
        elevator = variable(nc, "elem.joint.1004.Phi")[:, 1]
        fz = variable(nc, "elem.joint.23.F")[:, 2]
        my = variable(nc, "elem.joint.23.M")[:, 1]
    count = min(map(len, (time_values, pitch, elevator, fz, my)))
    mask = time_values[:count] >= SETTLED_START
    if np.count_nonzero(mask) < 10:
        raise RuntimeError(f"Caso incompleto: {perturbation.name}")
    return {
        "case": perturbation.name,
        "pitch_command_deg": perturbation.pitch_deg,
        "elevator_command_deg": perturbation.elevator_deg,
        "pitch_mean_deg": float(np.nanmean(pitch[:count][mask]) * RAD_TO_DEG),
        "elevator_mean_deg": float(np.nanmean(elevator[:count][mask]) * RAD_TO_DEG),
        "Fz_mean_lbf": float(np.nanmean(fz[:count][mask])),
        "Fz_std_lbf": float(np.nanstd(fz[:count][mask])),
        "My_mean_lbfin": float(np.nanmean(my[:count][mask])),
        "My_std_lbfin": float(np.nanstd(my[:count][mask])),
    }


def central_column(
    positive: dict[str, float | str],
    negative: dict[str, float | str],
    coordinate: str,
) -> np.ndarray:
    delta = math.radians(float(positive[coordinate]) - float(negative[coordinate]))
    if abs(delta) < 1.0e-12:
        raise RuntimeError(f"Perturbazione nulla per {coordinate}")
    return np.array(
        [
            (float(positive["Fz_mean_lbf"]) - float(negative["Fz_mean_lbf"])) / delta,
            (float(positive["My_mean_lbfin"]) - float(negative["My_mean_lbfin"])) / delta,
        ]
    )


def run() -> int:
    args = arguments()
    if args.pitch_step <= 0.0 or args.elevator_step <= 0.0 or args.jobs <= 0:
        raise ValueError("Passi e numero di job devono essere positivi")
    input_path = args.input.expanduser().resolve()
    output_directory = args.output_dir.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    template = input_path.read_text(encoding="utf-8")
    mbdyn = locate_mbdyn(args.mbdyn)

    scale = (1.146e-7 / REFERENCE_DENSITY) * (63.0 / REFERENCE_SPEED) ** 2
    pitch_center = 0.95525002 * scale - 0.80393236
    elevator_center = -0.19259962 * scale + 0.98569821
    perturbations = [
        Perturbation("baseline", pitch_center, elevator_center),
        Perturbation("pitch_plus", pitch_center + args.pitch_step, elevator_center),
        Perturbation("pitch_minus", pitch_center - args.pitch_step, elevator_center),
        Perturbation("elevator_plus", pitch_center, elevator_center + args.elevator_step),
        Perturbation("elevator_minus", pitch_center, elevator_center - args.elevator_step),
    ]

    print(f"Reference: V={REFERENCE_SPEED:g} m/s, rho={REFERENCE_DENSITY:.8e} IPS")
    print(f"Center: pitch={pitch_center:.6f} deg, BFL/BFR={elevator_center:.6f} deg")
    print(f"Output: {output_directory}")

    # ThreadPool is appropriate because each task waits on an external MBDyn process.
    from concurrent.futures import ThreadPoolExecutor, as_completed

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                launch_case,
                mbdyn,
                template,
                output_directory,
                perturbation,
                args.force,
            ): perturbation
            for perturbation in perturbations
        }
        for future in as_completed(futures):
            name, returncode, elapsed = future.result()
            print(f"{name}: returncode={returncode}, elapsed={elapsed:.1f} s", flush=True)
            if returncode != 0:
                raise RuntimeError(f"Calibrazione fallita: {name}")

    rows = [read_result(output_directory, item) for item in perturbations]
    by_name = {str(row["case"]): row for row in rows}
    theta_column = central_column(
        by_name["pitch_plus"], by_name["pitch_minus"], "pitch_mean_deg"
    )
    elevator_column = central_column(
        by_name["elevator_plus"],
        by_name["elevator_minus"],
        "elevator_mean_deg",
    )
    jacobian = np.column_stack((theta_column, elevator_column))
    inverse = np.linalg.inv(jacobian)
    condition_number = float(np.linalg.cond(jacobian))
    baseline = by_name["baseline"]
    residual = np.array(
        [float(baseline["Fz_mean_lbf"]), float(baseline["My_mean_lbfin"])]
    )
    center = np.radians(np.array([pitch_center, elevator_center]))
    newton_trim_deg = (center - inverse @ residual) * RAD_TO_DEG

    validation_case = Perturbation(
        "newton_validation", float(newton_trim_deg[0]), float(newton_trim_deg[1])
    )
    name, returncode, elapsed = launch_case(
        mbdyn, template, output_directory, validation_case, args.force
    )
    print(f"{name}: returncode={returncode}, elapsed={elapsed:.1f} s", flush=True)
    if returncode != 0:
        raise RuntimeError("Validazione Newton fallita")
    validation = read_result(output_directory, validation_case)
    rows.append(validation)

    with (output_directory / "calibration_cases.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "reference_speed_mps": REFERENCE_SPEED,
        "reference_density_ips": REFERENCE_DENSITY,
        "pitch_step_deg": args.pitch_step,
        "elevator_step_deg": args.elevator_step,
        "jacobian_rows": ["Fz_lbf", "My_lbfin"],
        "jacobian_columns": ["pitch_rad", "BFL_BFR_rad"],
        "jacobian": jacobian.tolist(),
        "inverse_jacobian": inverse.tolist(),
        "condition_number": condition_number,
        "baseline_residual": residual.tolist(),
        "one_step_newton_trim_deg": newton_trim_deg.tolist(),
        "validation_residual": [
            float(validation["Fz_mean_lbf"]),
            float(validation["My_mean_lbfin"]),
        ],
        "validation_standard_deviation": [
            float(validation["Fz_std_lbf"]),
            float(validation["My_std_lbfin"]),
        ],
    }
    (output_directory / "jacobian_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    print("\nJ = d[Fz, My]/d[pitch, BFL/BFR]")
    print(jacobian)
    print("\ninv(J) =")
    print(inverse)
    print(f"\ncond(J) = {condition_number:.6g}")
    print(
        "One-step Newton trim [pitch_deg, BFL/BFR_deg] = "
        f"[{newton_trim_deg[0]:.8f}, {newton_trim_deg[1]:.8f}]"
    )
    print(
        "Validation residual [Fz_lbf, My_lbfin] = "
        f"[{float(validation['Fz_mean_lbf']):.8f}, "
        f"{float(validation['My_mean_lbfin']):.8f}]"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
