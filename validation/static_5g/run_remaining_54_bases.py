#!/usr/bin/env python3
"""Resume the 5g MBDyn sweep for every elastic basis N=1,...,54."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

import netCDF4


ALL_COUNTS = tuple(range(1, 55))
REQUIRED_NODES = (990020, 991020)


def fem_mode_count(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if fields and fields[0].upper().startswith("REV") and len(fields) >= 3:
                return int(fields[2])
    raise RuntimeError(f"Cannot read FEM header: {path}")


def valid_result(path: Path) -> tuple[bool, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return False, "missing"
    try:
        with netCDF4.Dataset(path) as dataset:
            if "time" not in dataset.variables:
                return False, "time variable missing"
            time = dataset.variables["time"]
            if len(time) < 2 or float(time[-1]) < 5.9:
                return False, "simulation did not reach final time"
            for node in REQUIRED_NODES:
                key = f"node.struct.{node}.X"
                if key not in dataset.variables or dataset.variables[key].shape[0] != len(time):
                    return False, f"incomplete variable {key}"
    except (OSError, RuntimeError, ValueError) as error:
        return False, f"unreadable NetCDF: {error}"
    return True, "complete"


def run_case(mbdyn: str, case: Path, results: Path, environment: dict[str, str]) -> None:
    stem = case.stem
    output = results / stem
    run_log = results / f"{stem}.run.log"
    print(f"Running {case.name}", flush=True)
    with run_log.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(
            [mbdyn, "-f", case.name, "-o", str(output)],
            cwd=case.parent,
            stdout=stream,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"{case.name} failed; inspect {run_log}")
    valid, reason = valid_result(output.with_suffix(".nc"))
    if not valid:
        raise RuntimeError(f"{case.name} ended without a valid result: {reason}")


def main() -> int:
    study = Path(__file__).resolve().parent
    mbdyn_dir = study / "mbdyn"
    cases_dir = mbdyn_dir / "cases"
    results_dir = mbdyn_dir / "results"
    fem = mbdyn_dir / "mbdyn_modal.fem"

    parser = argparse.ArgumentParser(
        description="Run only missing MBDyn cases for all 54 elastic modal bases."
    )
    parser.add_argument("--mbdyn", default="/usr/local/mbdyn/bin/mbdyn")
    parser.add_argument(
        "--rerun", default="",
        help="Comma-separated mode counts to rerun even when their NetCDF is valid",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-analysis", action="store_true")
    args = parser.parse_args()

    if fem_mode_count(fem) != 60:
        raise RuntimeError(f"Expected a 60-mode FEM (6 rigid + 54 elastic): {fem}")
    rerun = {
        int(value) for value in re.split(r"[, ]+", args.rerun.strip()) if value
    }
    invalid_rerun = rerun - set(ALL_COUNTS)
    if invalid_rerun:
        raise ValueError(f"Invalid --rerun counts: {sorted(invalid_rerun)}")

    results_dir.mkdir(parents=True, exist_ok=True)
    status: dict[int, tuple[bool, str]] = {}
    for count in ALL_COUNTS:
        result = results_dir / f"gravity_5g_{count:02d}_elastic_modes.nc"
        status[count] = valid_result(result)
    completed = [count for count, (valid, _) in status.items() if valid and count not in rerun]
    pending = [count for count in ALL_COUNTS if count not in completed]

    print(f"Valid existing bases ({len(completed)}): {completed}")
    print(f"Bases still to run ({len(pending)}): {pending}")
    for count in pending:
        if status[count][1] != "missing":
            print(f"  N={count}: rerun required ({status[count][1]})")
    if args.dry_run:
        print("Dry run complete: no case generated or executed.")
        return 0

    subprocess.run(
        [
            sys.executable,
            str(mbdyn_dir / "generate_cases.py"),
            "--counts",
            ",".join(map(str, ALL_COUNTS)),
        ],
        cwd=mbdyn_dir,
        check=True,
    )

    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    for index, count in enumerate(pending, start=1):
        case = cases_dir / f"gravity_5g_{count:02d}_elastic_modes.mbd"
        if not case.is_file():
            raise RuntimeError(f"Generated case missing: {case}")
        print(f"[{index}/{len(pending)}]", end=" ", flush=True)
        run_case(args.mbdyn, case, results_dir, environment)

    final_invalid = []
    for count in ALL_COUNTS:
        valid, reason = valid_result(
            results_dir / f"gravity_5g_{count:02d}_elastic_modes.nc"
        )
        if not valid:
            final_invalid.append((count, reason))
    if final_invalid:
        raise RuntimeError(f"Sweep incomplete after execution: {final_invalid}")

    if not args.no_analysis:
        analysis_environment = environment.copy()
        analysis_environment.setdefault(
            "MPLCONFIGDIR", "/tmp/matplotlib-gravity-convergence"
        )
        subprocess.run(
            [sys.executable, str(study / "analyze_convergence.py")],
            cwd=study,
            env=analysis_environment,
            check=True,
        )
    print("All 54 elastic bases are complete; convergence outputs regenerated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
