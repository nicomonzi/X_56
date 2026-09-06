#!/usr/bin/env python3
"""Generate and run every MBDyn gravity convergence case."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--mbdyn", default="/usr/local/mbdyn/bin/mbdyn")
    parser.add_argument("--counts", default=None, help="Optional comma-separated elastic mode counts")
    parser.add_argument("--generate-only", action="store_true")
    args = parser.parse_args()

    generate = [sys.executable, str(root / "generate_cases.py")]
    if args.counts:
        generate += ["--counts", args.counts]
    subprocess.run(generate, cwd=root, check=True)
    if args.generate_only:
        return 0

    result_dir = root / "results"
    result_dir.mkdir(exist_ok=True)
    for pattern in (
        "gravity_5g_*_elastic_modes.nc",
        "gravity_5g_*_elastic_modes.log",
        "gravity_5g_*_elastic_modes.out",
        "gravity_5g_*_elastic_modes.run.log",
    ):
        for stale_result in result_dir.glob(pattern):
            stale_result.unlink()
    cases = sorted((root / "cases").glob("gravity_5g_*_elastic_modes.mbd"))
    for index, case in enumerate(cases, start=1):
        stem = case.stem
        output = result_dir / stem
        log = result_dir / f"{stem}.run.log"
        print(f"[{index}/{len(cases)}] {case.name}", flush=True)
        with log.open("w", encoding="utf-8") as stream:
            subprocess.run(
                [args.mbdyn, "-f", case.name, "-o", str(output)],
                cwd=case.parent,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=True,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
