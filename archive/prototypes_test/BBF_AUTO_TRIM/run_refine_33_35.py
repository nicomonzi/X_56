#!/usr/bin/env python3
"""Run six targeted cases around the nonlinear 33 and 35 m/s results."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parent
MAIN = ROOT / "main_bbf.mbd"
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")
OUTPUT = Path("/mnt/c/Users/Utente/Desktop/RESULTS")

# Sei soli punti intermedi: tre tra 33--34 e tre tra 34--35 m/s.
# La risoluzione di 0.25 m/s e' sufficiente per localizzare le due transizioni
# senza lanciare la griglia fitta e le bisezioni della prima versione.
REFINEMENT_VELOCITIES = (33.25, 33.50, 33.75, 34.25, 34.50, 34.75)

INCH_TO_M = 0.0254
RAD_TO_DEG = 180.0 / math.pi
BASE_NODE = 990000
TRIM_PITCH_PID = 9101
TRIM_END = 12.0
BFF_START = 13.75
BFF_END = 28.75


def velocity_key(velocity: float) -> int:
    """Integer micro-m/s key: avoids duplicate runs from float roundoff."""
    return int(round(velocity * 1_000_000))


def case_directory(velocity: float) -> Path:
    """Use decimal-safe names distinct from the existing integer cases."""
    if abs(velocity - round(velocity)) < 1.0e-9:
        return OUTPUT / f"V_{int(round(velocity)):03d}_mps"
    token = f"{velocity:07.3f}".replace(".", "p")
    return OUTPUT / f"V_{token}_mps"


def nc_data(dataset: Dataset, name: str) -> np.ndarray:
    if name not in dataset.variables:
        raise KeyError(f"Missing NetCDF variable {name}")
    return np.ma.filled(dataset.variables[name][:], np.nan).astype(float)


def classify(result: Path) -> dict[str, float | bool | str]:
    """Apply the same rigid-body validity limits used by analyze_sweep.py."""
    with Dataset(result) as dataset:
        time = nc_data(dataset, "time")
        if len(time) == 0 or time[-1] < BFF_END:
            raise RuntimeError(f"Incomplete result, preserved without changes: {result}")
        position = nc_data(dataset, f"node.struct.{BASE_NODE}.X")
        attitude = nc_data(dataset, f"node.struct.{BASE_NODE}.Phi")
        trim_pitch = (
            nc_data(dataset, f"elem.loadable.{TRIM_PITCH_PID}.output").squeeze()
            * RAD_TO_DEG
        )

    trim_window = (time >= TRIM_END - 1.0) & (time < TRIM_END)
    bff_window = (time >= BFF_START + 0.02) & (time <= BFF_END)
    if np.count_nonzero(trim_window) < 10 or np.count_nonzero(bff_window) < 100:
        raise RuntimeError(f"Required trim/BFF windows are absent: {result}")

    heave = (position[:, 2] - position[0, 2]) * INCH_TO_M
    pitch = attitude[:, 1] * RAD_TO_DEG
    trim_pitch_mean = float(np.mean(trim_pitch[trim_window]))
    heave_range = float(np.ptp(heave[bff_window]))
    pitch_deviation = float(
        np.max(np.abs(pitch[bff_window] - trim_pitch_mean))
    )
    divergent = bool(heave_range > 5.0 or pitch_deviation > 10.0)
    return {
        "divergent": divergent,
        "state": "DIVERGENT" if divergent else "BOUNDED",
        "heave_range_m": heave_range,
        "pitch_deviation_deg": pitch_deviation,
    }


def create_input(velocity: float, case: Path) -> tuple[Path, str]:
    source = MAIN.read_text(encoding="utf-8")
    pattern = r"(?m)^set:\s*const\s+real\s+V_INF\s*=\s*[^;]+;"
    case_source, replacements = re.subn(
        pattern,
        f"set: const real V_INF = {velocity:.6f};",
        source,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("Could not replace the unique V_INF definition")

    case.mkdir(parents=True, exist_ok=True)
    permanent_input = case / "case_input.mbd"
    if permanent_input.exists():
        old_source = permanent_input.read_text(encoding="utf-8")
        if old_source != case_source:
            raise RuntimeError(
                f"Existing input differs and will not be overwritten: {permanent_input}"
            )
    else:
        permanent_input.write_text(case_source, encoding="utf-8")
    return permanent_input, case_source


def run_or_reuse(velocity: float) -> dict[str, float | bool | str]:
    case = case_directory(velocity)
    result = case / "case.nc"

    if result.exists():
        metrics = classify(result)
        print(
            f"Reuse {velocity:.6f} m/s: {metrics['state']}, "
            f"heave={metrics['heave_range_m']:.4g} m, "
            f"pitch={metrics['pitch_deviation_deg']:.4g} deg",
            flush=True,
        )
        return metrics

    if case.exists() and any(case.iterdir()):
        raise RuntimeError(
            f"Non-empty case without case.nc will not be overwritten: {case}"
        )

    _, case_source = create_input(velocity, case)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=".run_refine_",
        suffix=".mbd",
        dir=ROOT,
        delete=False,
    ) as temporary:
        temporary.write(case_source)
        active_input = Path(temporary.name)

    print(f"\n=== Refinement run: V_INF = {velocity:.6f} m/s ===", flush=True)
    try:
        subprocess.run(
            [
                str(MBDYN),
                "-f",
                str(active_input),
                "-o",
                str(case / "case"),
            ],
            cwd=ROOT,
            check=True,
        )
    finally:
        active_input.unlink(missing_ok=True)

    metrics = classify(result)
    print(
        f"Result {velocity:.6f} m/s: {metrics['state']}, "
        f"heave={metrics['heave_range_m']:.4g} m, "
        f"pitch={metrics['pitch_deviation_deg']:.4g} deg",
        flush=True,
    )
    return metrics


def write_summary(results: dict[int, tuple[float, dict]]) -> None:
    destination = OUTPUT / "refinement_33_35_summary.csv"
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "velocity_mps",
                "state",
                "heave_range_m",
                "pitch_deviation_deg",
                "result_directory",
            ]
        )
        for velocity, metrics in sorted(results.values()):
            writer.writerow(
                [
                    f"{velocity:.6f}",
                    metrics["state"],
                    metrics["heave_range_m"],
                    metrics["pitch_deviation_deg"],
                    case_directory(velocity),
                ]
            )
    print(f"\nRefinement summary: {destination}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run six targeted BBF_AUTO_TRIM cases around the existing "
            "nonlinear results at 33 and 35 m/s."
        )
    )
    parser.parse_args()

    if not MAIN.is_file() or not MBDYN.is_file():
        raise SystemExit("main_bbf.mbd or the MBDyn executable is missing")

    results: dict[int, tuple[float, dict]] = {}

    def evaluate(velocity: float) -> dict:
        key = velocity_key(velocity)
        if key not in results:
            results[key] = (velocity, run_or_reuse(velocity))
            write_summary(results)
        return results[key][1]

    # Existing integer cases are read, never regenerated or overwritten.
    evaluate(33.0)
    evaluate(34.0)
    evaluate(35.0)

    # Exactly six new runs, all written below the RESULTS directory on C:.
    for velocity in REFINEMENT_VELOCITIES:
        evaluate(velocity)

    ordered = sorted(results.values())
    transitions: list[tuple[float, float]] = []
    for (left_v, left_m), (right_v, right_m) in zip(ordered, ordered[1:]):
        if left_m["divergent"] != right_m["divergent"]:
            transitions.append((left_v, right_v))

    if not transitions:
        print(
            "\nNo bounded/divergent change was found on the 0.25 m/s grid "
            "between 33 and 35 m/s."
        )
        return

    print("\nBounded/divergent transition brackets:", flush=True)
    for left, right in transitions:
        left_state = results[velocity_key(left)][1]["state"]
        right_state = results[velocity_key(right)][1]["state"]
        print(
            f"  {left:.6f} ({left_state})  --  "
            f"{right:.6f} m/s ({right_state})",
            flush=True,
        )
    print(
        "\nThese 0.25 m/s brackets identify rigid-body departure boundaries, not by "
        "themselves a certified flutter speed. Re-run analyze_sweep.py to "
        "check modal damping and fit quality.",
        flush=True,
    )


if __name__ == "__main__":
    main()
