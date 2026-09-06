#!/usr/bin/env python3
"""End-to-end 5g workflow: femgen, MBDyn modal sweep, and final analysis."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def fem_header(path: Path) -> tuple[int, int]:
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if fields and fields[0].upper().startswith("REV") and len(fields) >= 3:
                return int(fields[1]), int(fields[2])
    raise RuntimeError(f"Cannot read FEM header: {path}")


def validate_fem(path: Path, expected_modes: int = 60) -> tuple[int, int]:
    nodes, modes = fem_header(path)
    text_tail = ""
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "RECORD GROUP 10" in line or "RECORD GROUP 11" in line:
                text_tail += line
    if modes != expected_modes:
        raise RuntimeError(f"Expected {expected_modes} FEM modes, found {modes}")
    if "RECORD GROUP 10" not in text_tail or "RECORD GROUP 11" not in text_tail:
        raise RuntimeError("Incomplete FEM: modal stiffness or lumped mass record missing")
    return nodes, modes


def run_case(
    mbdyn: str, case: Path, result_dir: Path, environment: dict[str, str]
) -> tuple[str, float]:
    stem = case.stem
    output = result_dir / stem
    log = result_dir / f"{stem}.run.log"
    command = [mbdyn, "-f", case.name, "-o", str(output)]
    with log.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(
            command,
            cwd=case.parent,
            stdout=stream,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"{case.name} failed; inspect {log}")
    mode_match = re.search(r"_(\d+)_elastic_modes$", stem)
    count = float(mode_match.group(1)) if mode_match else float("nan")
    return case.name, count


def main() -> int:
    study = Path(__file__).resolve().parent
    package = study.parent
    modal_main = package / "01_SOL103_60_MODES/MAIN"
    mbdyn_dir = study / "mbdyn"
    result_dir = mbdyn_dir / "results"

    parser = argparse.ArgumentParser(
        description="Generate the 60-mode FEM, run the 5g MBDyn sweep, and analyze it."
    )
    parser.add_argument("--femgen", default="/usr/local/mbdyn/bin/femgen")
    parser.add_argument("--mbdyn", default="/usr/local/mbdyn/bin/mbdyn")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--counts", default=None, help="Comma-separated elastic mode counts")
    parser.add_argument("--skip-femgen", action="store_true")
    parser.add_argument("--absolute-threshold", type=float, default=0.01)
    parser.add_argument("--relative-threshold", type=float, default=0.01)
    args = parser.parse_args()
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")

    generated_fem = modal_main / "mbdyn_modal_60.fem"
    installed_fem = mbdyn_dir / "mbdyn_modal.fem"
    if not args.skip_femgen:
        print("[1/4] Running femgen on the corrected 60-mode OP2/MAT", flush=True)
        subprocess.run(
            [args.femgen, "sol103_60_modes", "-o", "mbdyn_modal_60"],
            cwd=modal_main,
            check=True,
        )
        nodes, modes = validate_fem(generated_fem)
        shutil.copy2(generated_fem, installed_fem)
        print(f"      Installed FEM: {nodes} nodes, {modes} total modes", flush=True)
    else:
        nodes, modes = validate_fem(installed_fem)
        print(f"[1/4] Using installed FEM: {nodes} nodes, {modes} total modes", flush=True)

    print("[2/4] Generating modal-truncation cases", flush=True)
    generate = [sys.executable, str(mbdyn_dir / "generate_cases.py")]
    if args.counts:
        generate += ["--counts", args.counts]
    subprocess.run(generate, cwd=mbdyn_dir, check=True)

    result_dir.mkdir(parents=True, exist_ok=True)
    for pattern in (
        "gravity_5g_*_elastic_modes.nc",
        "gravity_5g_*_elastic_modes.log",
        "gravity_5g_*_elastic_modes.out",
        "gravity_5g_*_elastic_modes.run.log",
    ):
        for stale in result_dir.glob(pattern):
            stale.unlink()
    cases = sorted((mbdyn_dir / "cases").glob("gravity_5g_*_elastic_modes.mbd"))
    if not cases:
        raise RuntimeError("No MBDyn cases were generated")

    print(f"[3/4] Running {len(cases)} MBDyn cases with {args.jobs} parallel jobs", flush=True)
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    completed_count = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_case, args.mbdyn, case, result_dir, environment): case
            for case in cases
        }
        for future in as_completed(futures):
            name, _ = future.result()
            completed_count += 1
            print(f"      [{completed_count}/{len(cases)}] completed {name}", flush=True)

    print("[4/4] Comparing MBDyn tips with the Nastran SOL 101 reference", flush=True)
    analysis = [
        sys.executable,
        str(study / "analyze_convergence.py"),
        "--absolute-threshold",
        str(args.absolute_threshold),
        "--relative-threshold",
        str(args.relative_threshold),
    ]
    analysis_environment = environment.copy()
    analysis_environment.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-gravity-convergence")
    subprocess.run(analysis, cwd=study, env=analysis_environment, check=True)
    print(f"Workflow complete. Plots and table: {study / 'plots'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
