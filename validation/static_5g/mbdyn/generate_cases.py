#!/usr/bin/env python3
"""Generate MBDyn 5g cases using progressively richer elastic modal bases."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


DEFAULT_COUNTS = (1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 30, 40, 50, 54)


def fem_mode_count(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if line.lstrip().startswith("**") or not line.strip():
                continue
            fields = line.split()
            if fields[0].upper().startswith("REV") and len(fields) >= 3:
                return int(fields[2])
    raise ValueError(f"Cannot read modal count from {path}")


def parse_counts(specification: str) -> list[int]:
    values = sorted({int(value) for value in re.split(r"[, ]+", specification) if value})
    if not values or values[0] < 1:
        raise ValueError("Mode counts must be positive integers")
    return values


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--fem", type=Path, default=root / "mbdyn_modal.fem")
    parser.add_argument("--template", type=Path, default=root / "gravity_5g_template.mbd")
    parser.add_argument("--output", type=Path, default=root / "cases")
    parser.add_argument("--counts", default=",".join(map(str, DEFAULT_COUNTS)))
    args = parser.parse_args()

    total_modes = fem_mode_count(args.fem)
    elastic_modes = total_modes - 6
    if elastic_modes < 1:
        raise ValueError("The FEM must contain six rigid-body modes followed by elastic modes")
    requested = parse_counts(args.counts)
    counts = [count for count in requested if count <= elastic_modes]
    if not counts:
        raise ValueError(f"No requested count fits the {elastic_modes} elastic FEM modes")
    if counts[-1] != elastic_modes:
        counts.append(elastic_modes)

    template = args.template.read_text(encoding="utf-8")
    args.output.mkdir(parents=True, exist_ok=True)
    for old_case in args.output.glob("gravity_5g_*_elastic_modes.mbd"):
        old_case.unlink()
    for count in counts:
        mode_list = ", ".join(str(mode) for mode in range(7, 7 + count))
        content = template.replace("@MODE_COUNT@", str(count)).replace(
            "@MODE_LIST@", mode_list
        )
        path = args.output / f"gravity_5g_{count:02d}_elastic_modes.mbd"
        path.write_text(content, encoding="utf-8")
        print(path.name)
    print(f"Generated {len(counts)} cases from a {total_modes}-mode FEM ({elastic_modes} elastic).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
