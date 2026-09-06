#!/usr/bin/env python3
"""Run coupled MBDyn/DUST velocities sequentially; no analysis is started."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=float, default=23.0, help="first m/s")
    parser.add_argument("--stop", type=float, default=36.0, help="last m/s")
    parser.add_argument("--step", type=float, default=1.0, help="velocity step")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    if args.step <= 0.0 or args.stop < args.start:
        raise SystemExit("Require step > 0 and stop >= start")

    count = int(round((args.stop - args.start) / args.step))
    velocities = [args.start + index * args.step for index in range(count + 1)]
    for velocity in velocities:
        print(f"\n=== Coupled case: {velocity:g} m/s ===", flush=True)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "run_coupled.py"),
                "--velocity",
                f"{velocity:.12g}",
            ],
            cwd=ROOT,
            check=True,
        )


if __name__ == "__main__":
    main()
