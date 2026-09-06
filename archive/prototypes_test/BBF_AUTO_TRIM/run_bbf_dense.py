#!/usr/bin/env python3
"""Run a dense, localized sweep around the suspected BBF boundary."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "run_bbf_boundary.py"
DEFAULT_OUTPUT = Path(
    "/mnt/c/Users/Utente/Desktop/RESULTS_BBF_DT002_MODES_7_23_DENSE"
)

# MODIFICA: lo sweep precedente colloca lo zero del ramo pitch--FEM7 vicino
# a 42.38 m/s. Si riusano i casi esistenti a 42 e 43 m/s e si campiona ogni
# 0.05 m/s nella zona dello zero, con pochi punti esterni di controllo.
VELOCITIES = (
    42.00,
    42.10,
    42.20,
    42.25,
    42.30,
    42.35,
    42.40,
    42.45,
    42.50,
    42.55,
    42.60,
    42.70,
    42.80,
    42.90,
    43.00,
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run 15 BBF cases around the estimated 42.38 m/s boundary, "
            "using 0.05 m/s resolution near the damping zero."
        )
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Result root; default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--no-analysis",
        action="store_true",
        help="Run MBDyn only and skip the final V-g/V-f analysis.",
    )
    return parser.parse_args()


def main() -> None:
    args = arguments()
    command = [
        sys.executable,
        str(RUNNER),
        "--output",
        str(args.output.expanduser().resolve()),
        "--velocities",
        *(f"{velocity:.2f}" for velocity in VELOCITIES),
    ]
    if args.no_analysis:
        command.append("--no-analysis")
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
