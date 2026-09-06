#!/usr/bin/env python3
"""Build the effective X-56 C81 used by the strip-aerodynamic model.

The correction is deliberately limited to quantities that can be supported by
the available full-aircraft data: effective lift slope and a small section
pitch-moment slope.  Profile drag from the original NACA 0012 table is retained.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "INCLUDE/naca0012.c81"
DEFAULT_OUTPUT = ROOT / "INCLUDE/x56_effective.c81"
TABLE_SHAPES = ((11, 39), (11, 65), (9, 47))


@dataclass
class Table:
    mach: list[float]
    alpha: list[float]
    values: list[list[float]]


def fields(line: str) -> list[str]:
    padded = line.rstrip("\n").ljust(80)
    return [padded[i : i + 7].strip() for i in range(0, 77, 7)]


def read_c81(path: Path) -> tuple[str, list[Table]]:
    lines = path.read_text().splitlines()
    header = lines[0]
    cursor = 1
    tables = []
    for n_mach, n_alpha in TABLE_SHAPES:
        mach: list[float] = []
        while len(mach) < n_mach:
            row = fields(lines[cursor]); cursor += 1
            mach.extend(float(value) for value in row[1:] if value)
        alpha: list[float] = []
        values: list[list[float]] = []
        for _ in range(n_alpha):
            row = fields(lines[cursor]); cursor += 1
            if not row[0]:
                raise ValueError(f"missing alpha at line {cursor} in {path}")
            angle = float(row[0])
            coefficients = [float(value) for value in row[1:] if value]
            while len(coefficients) < n_mach:
                row = fields(lines[cursor]); cursor += 1
                if row[0]:
                    raise ValueError(f"unexpected new alpha at line {cursor} in {path}")
                coefficients.extend(float(value) for value in row[1:] if value)
            alpha.append(angle); values.append(coefficients[:n_mach])
        tables.append(Table(mach[:n_mach], alpha, values))
    return header, tables


def number(value: float) -> str:
    if abs(value) < 5e-12:
        value = 0.0
    result = f"{value:.5g}"
    if result == "0":
        return "0.".ljust(7)
    if "." not in result and "e" not in result.lower():
        result += "."
    if result.startswith("0."):
        result = result[1:]
    elif result.startswith("-0."):
        result = "-" + result[2:]
    if len(result) > 7:
        result = f"{value:.1e}"
    if len(result) > 7:
        raise ValueError(f"C81 value does not fit a seven-character field: {value}")
    return result.ljust(7)


def line(values: list[float], first: float | None = None) -> str:
    prefix = "       " if first is None else number(first)
    return (prefix + "".join(number(value) for value in values)).ljust(80)


def write_c81(path: Path, header: str, tables: list[Table]) -> None:
    output = [header[:80].ljust(80)]
    for table in tables:
        output.append(line(table.mach[:9]))
        if len(table.mach) > 9:
            output.append(line(table.mach[9:]))
        for angle, values in zip(table.alpha, table.values):
            output.append(line(values[:9], angle))
            if len(values) > 9:
                output.append(line(values[9:]))
    path.write_text("\n".join(output) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cl-scale", type=float, default=0.773200)
    parser.add_argument("--cm-alpha", type=float, default=-0.005930,
                        help="added local Cm slope per degree for |alpha|<=16 deg")
    args = parser.parse_args()
    header, tables = read_c81(args.source)
    lift, _drag, moment = tables
    lift.values = [[args.cl_scale * value for value in row] for row in lift.values]
    for i, angle in enumerate(moment.alpha):
        magnitude = abs(angle)
        taper = 1.0 if magnitude <= 16.0 else max(0.0, (30.0 - magnitude) / 14.0)
        moment.values[i] = [value + args.cm_alpha * angle * taper for value in moment.values[i]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    calibrated_header = "X-56 EFFECTIVE C81".ljust(28) + header[28:]
    write_c81(args.output, calibrated_header, tables)
    # Parse the result again: catches malformed fixed-width output before MBDyn.
    read_c81(args.output)
    print(f"wrote {args.output}: CL scale={args.cl_scale:.6f}, Cm_alpha={args.cm_alpha:+.6f}/deg")


if __name__ == "__main__":
    main()
