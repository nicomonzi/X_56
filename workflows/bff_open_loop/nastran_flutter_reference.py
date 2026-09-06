#!/usr/bin/env python3
"""Extract the released X-56A SOL 145 flutter branch used as BFF reference."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

KNOT_TO_MPS = 0.5144444444444445
INCH_TO_METER = 0.0254
VREF_INPS_PER_KEAS = 21.9849
REFINED_F06 = Path(
    "/home/nicomonzi/ZENO/X56_NASTRAN_BFF_REFINEMENT/FLUTTER_TEST/"
    "x56_bff_refined.f06"
)
REPO_ROOT = Path(__file__).resolve().parents[2]
RELEASED_F06 = REPO_ROOT / "models/nastran/x56_r11/FLUTTER_TEST/nsfluttr_cfg24611_10lb.f06"
DEFAULT_F06 = REFINED_F06 if REFINED_F06.exists() else RELEASED_F06


def parse_flutter_point(path: Path = DEFAULT_F06, point: int = 7) -> list[dict[str, float]]:
    """Read one PK flutter-summary branch from the released F06 file.

    The deck uses the density-corrected ``PARAM,VREF,21.9849`` to print
    velocity in KEAS; the FLFACT values themselves are TAS in inch/s.  MBDyn
    uses TAS, so ``velocity_mps`` below deliberately means TAS.  The complex-eigenvalue real
    part is sigma directly and is preferable to reconstructing it from the
    rounded NASTRAN damping column.
    """
    lines = path.read_text(errors="ignore").splitlines()
    header = re.compile(rf"\bPOINT\s*=\s*{point}\b")
    start = next((i for i, line in enumerate(lines) if header.search(line)), None)
    if start is None:
        raise ValueError(f"POINT={point} not found in {path}")
    number = r"[-+]?\d*\.\d+(?:[ED][-+]?\d+)?"
    row_re = re.compile(rf"^\s*({number})(?:\s+({number})){{6}}\s*$")
    rows: list[dict[str, float]] = []
    for line in lines[start + 1:]:
        if rows and "FLUTTER  SUMMARY" in line:
            break
        tokens = line.replace("D", "E").split()
        if len(tokens) != 7 or not row_re.match(line.replace("D", "E")):
            continue
        kfreq, inv_kfreq, velocity_keas, damping, frequency, sigma, omega = map(float, tokens)
        rows.append({
            "point": float(point),
            "velocity_keas": velocity_keas,
            "velocity_eas_mps": velocity_keas * KNOT_TO_MPS,
            "velocity_mps": velocity_keas * VREF_INPS_PER_KEAS * INCH_TO_METER,
            "frequency_hz": frequency,
            "sigma_per_s": sigma,
            "nastran_damping_g": damping,
            "kfreq": kfreq,
            "inverse_kfreq": inv_kfreq,
            "omega_rad_s": omega,
        })
    if not rows:
        raise ValueError(f"no flutter rows found for POINT={point} in {path}")
    return rows


def zero_crossing(rows: list[dict[str, float]]) -> dict[str, float] | None:
    for lower, upper in zip(rows, rows[1:]):
        a, b = lower["sigma_per_s"], upper["sigma_per_s"]
        if a <= 0.0 < b:
            fraction = -a / (b - a)
            return {
                "velocity_mps": lower["velocity_mps"] + fraction * (upper["velocity_mps"] - lower["velocity_mps"]),
                "velocity_keas": lower["velocity_keas"] + fraction * (upper["velocity_keas"] - lower["velocity_keas"]),
                "velocity_eas_mps": lower["velocity_eas_mps"] + fraction * (upper["velocity_eas_mps"] - lower["velocity_eas_mps"]),
                "frequency_hz": lower["frequency_hz"] + fraction * (upper["frequency_hz"] - lower["frequency_hz"]),
            }
    return None


def interpolate_at_velocity(rows: list[dict[str, float]], velocity_mps: float) -> dict[str, float] | None:
    """Linearly interpolate one branch without extrapolating it."""
    ordered = sorted(rows, key=lambda row: row["velocity_mps"])
    for lower, upper in zip(ordered, ordered[1:]):
        if lower["velocity_mps"] <= velocity_mps <= upper["velocity_mps"]:
            span = upper["velocity_mps"] - lower["velocity_mps"]
            ratio = (velocity_mps - lower["velocity_mps"]) / span
            return {
                "velocity_mps": float(velocity_mps),
                "frequency_hz": lower["frequency_hz"] + ratio * (upper["frequency_hz"] - lower["frequency_hz"]),
                "sigma_per_s": lower["sigma_per_s"] + ratio * (upper["sigma_per_s"] - lower["sigma_per_s"]),
                "lower_grid_mps": lower["velocity_mps"],
                "upper_grid_mps": upper["velocity_mps"],
            }
    return None


def write_reference_files(output: Path, path: Path = DEFAULT_F06, point: int = 7) -> dict:
    rows = parse_flutter_point(path, point)
    payload = {
        "source": str(path),
        "point": point,
        "interpretation": "BFF branch: rigid short-period coupled with first symmetric wing bending",
        "velocity_note": "velocity_mps is TAS; F06 VELOCITY is KEAS because VREF is density corrected",
        "zero_crossing_linear": zero_crossing(rows),
        "rows": rows,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "nastran_bff_reference.json").write_text(json.dumps(payload, indent=2) + "\n")
    with (output / "nastran_bff_reference.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--f06", type=Path, default=DEFAULT_F06)
    parser.add_argument("--point", type=int, default=7)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = parse_flutter_point(args.f06, args.point)
    payload = {
        "source": str(args.f06), "point": args.point,
        "interpretation": "BFF branch: rigid short-period coupled with first symmetric wing bending",
        "velocity_note": "velocity_mps is TAS; F06 VELOCITY is KEAS because VREF is density corrected",
        "zero_crossing_linear": zero_crossing(rows), "rows": rows,
    }
    if args.output:
        payload = write_reference_files(args.output, args.f06, args.point)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
