#!/usr/bin/env python3
"""NASA-style full-surface-hold SAS-off sweep from 50 to 70 m/s (manual launch only).

Outputs default to C:\\Users\\Utente\\Desktop\\BFF_open_loop and use the
NASA_OL_V_* prefix.  The output directory is expected to be empty at launch.
"""

from __future__ import annotations

import argparse
import json
import os
from decimal import Decimal
from pathlib import Path

from analyse_open_loop import analyze_case, write_global_plots
from nastran_flutter_reference import write_reference_files
from run_case import DEFAULT_OUTPUT, OUTPUT_ENV, run_case


def clean_generated_outputs(output: Path) -> None:
    """Remove only files produced by this sweep, never unrelated user files."""
    output.mkdir(parents=True, exist_ok=True)
    patterns = (
        "NASA_OL_V_*",
        "sweep_summary.json",
        "local_candidate_overview.png",
        "nastran_bff_point7.*",
    )
    removed = 0
    for pattern in patterns:
        for path in output.glob(pattern):
            if path.is_file():
                path.unlink()
                removed += 1
    plots = output / "plots"
    if plots.is_dir():
        for path in plots.glob("NASA_OL_V_*.png"):
            if path.is_file():
                path.unlink()
                removed += 1
    print(f"[clean] rimossi {removed} file generati da run precedenti")


def grid(start: float, stop: float, step: float) -> list[float]:
    a, b, h = Decimal(str(start)), Decimal(str(stop)), Decimal(str(step))
    if h <= 0 or b < a:
        raise ValueError("richiesti start <= stop e step > 0")
    values: list[float] = []
    while a <= b + Decimal("1e-12"):
        values.append(float(a))
        a += h
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=float, default=50.0)
    parser.add_argument("--stop", type=float, default=70.0)
    parser.add_argument("--step", type=float, default=2.5)
    parser.add_argument("--refine-tolerance", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=Path(os.environ.get(OUTPUT_ENV, DEFAULT_OUTPUT)))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--clean", action="store_true",
        help="remove previous sweep-generated files from the output directory before launch",
    )
    args = parser.parse_args()

    if args.clean:
        clean_generated_outputs(args.output)

    summaries: dict[float, dict] = {}

    def evaluate(v: float) -> dict:
        nc = run_case(v, args.output, args.overwrite)
        tracked_neighbors = [
            (speed, value) for speed, value in summaries.items()
            if value.get("matrix_pencil_candidate", {}).get("frequency_hz") is not None
        ]
        reference_frequency = None
        if tracked_neighbors:
            _, neighbor = min(tracked_neighbors, key=lambda item: abs(item[0] - v))
            reference_frequency = float(neighbor["matrix_pencil_candidate"]["frequency_hz"])
        summary = analyze_case(
            nc, args.output,
            tracking_reference_frequency_hz=reference_frequency,
        )
        if summary.get("identification_valid") and tracked_neighbors:
            neighbor_speed, neighbor = min(tracked_neighbors, key=lambda item: abs(item[0] - v))
            frequency = float(summary["frequency_swb1_hz"])
            neighbor_frequency = float(neighbor["frequency_swb1_hz"])
            tolerance = max(0.40, 0.15 * 0.5 * (frequency + neighbor_frequency))
            summary["tracking_reference_velocity_mps"] = neighbor_speed
            summary["tracking_frequency_difference_hz"] = abs(frequency - neighbor_frequency)
            summary["tracking_valid"] = abs(frequency - neighbor_frequency) <= tolerance
        else:
            summary["tracking_reference_velocity_mps"] = None
            summary["tracking_frequency_difference_hz"] = None
            summary["tracking_valid"] = bool(summary.get("identification_valid"))
        summaries[v] = summary
        return summary

    for velocity in grid(args.start, args.stop, args.step):
        evaluate(velocity)

    valid = sorted(
        (v, s["sigma_swb1_per_s"])
        for v, s in summaries.items()
        if s["identification_valid"] and s.get("tracking_valid", False)
        and s["sigma_swb1_per_s"] is not None
    )
    brackets = [(a, b) for a, b in zip(valid, valid[1:]) if a[1] <= 0.0 < b[1]]
    if brackets:
        lo, hi = brackets[0][0][0], brackets[0][1][0]
        while hi - lo > args.refine_tolerance:
            mid = 0.5 * (lo + hi)
            midpoint = evaluate(mid)
            sigma = midpoint["sigma_swb1_per_s"]
            if sigma is None or not midpoint.get("tracking_valid", False):
                break
            if sigma > 0.0:
                hi = mid
            else:
                lo = mid

    ordered = [summaries[v] for v in sorted(summaries)]
    (args.output / "sweep_summary.json").write_text(json.dumps(ordered, indent=2))
    write_reference_files(args.output)
    write_global_plots(ordered, args.output)


if __name__ == "__main__":
    main()
