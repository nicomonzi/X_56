#!/usr/bin/env python3
"""Prepare or explicitly run one paired MBDyn campaign."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from netCDF4 import Dataset

from campaign import BASELINE_DIR, CONFIG, MODEL, ROOT, build_cases, render_case, source_fingerprints

DEFAULT_OUTPUT_ROOT = Path(CONFIG["execution"]["output_root"])


def result_is_complete(nc: Path, mbd: Path) -> bool:
    """Reject a truncated NetCDF instead of silently treating it as cached."""
    if not nc.is_file():
        return False
    text = mbd.read_text()
    values = {}
    for name in ("SAS_OFF_START", "SAS_OFF_DURATION"):
        match = re.search(
            rf"(?m)^set:\s*const\s+real\s+{name}\s*=\s*([-+0-9.eE]+)\s*;",
            text,
        )
        if not match:
            return False
        values[name] = float(match.group(1))
    try:
        with Dataset(nc) as data:
            time = data["time"][:]
        return bool(
            len(time) >= 50
            and float(time[-1])
            >= values["SAS_OFF_START"] + values["SAS_OFF_DURATION"] - 0.011
        )
    except Exception:
        return False


def write_manifest(rows: list[dict], path: Path, fingerprints: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    path.with_suffix(".json").write_text(json.dumps({
        "source_sha256": fingerprints,
        "cases": rows,
    }, indent=2))


def prepare(
    campaign_name: str, output: Path, overwrite: bool = False,
) -> tuple[list, list[dict], Path]:
    cases = build_cases(campaign_name)
    case_dir = output / "cases"
    case_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for index, case in enumerate(cases, start=1):
        mbd = case_dir / f"{case.stem}.mbd"
        nc = mbd.with_suffix(".nc")
        rendered = render_case(case)
        if mbd.is_file() and nc.is_file() and mbd.read_text() != rendered and not overwrite:
            raise RuntimeError(
                f"input cambiato ma risultato esistente: {nc}; "
                "usare --overwrite per rigenerare ed eseguire coerentemente"
            )
        mbd.write_text(rendered)
        metadata = case.metadata()
        already_complete = result_is_complete(nc, mbd) and not overwrite
        rows.append({
            "case_index": index,
            "campaign": campaign_name,
            "velocity_mps": case.velocity_mps,
            "nominal_load_factor": case.nominal_load_factor,
            "time_step_s": case.time_step_s,
            "mode7_frequency_scale": case.mode7_frequency_scale,
            "excited": case.excited,
            "pitch_rate_command_deg_s": metadata["pitch_rate_command_deg_s"],
            "pitch_amplitude_deg": metadata["pitch_amplitude_deg"],
            "stiffness_delta_k7": metadata["stiffness_delta_k7"],
            "stem": case.stem,
            "input_sha256": hashlib.sha256(rendered.encode()).hexdigest(),
            "status": "existing" if already_complete else "prepared",
            "message": (
                "existing result with identical input"
                if already_complete else "input MBDyn generated; no solver run"
            ),
        })
    manifest = output / "manifest.csv"
    write_manifest(rows, manifest, source_fingerprints())
    return cases, rows, manifest


def run_one(case, case_dir: Path, overwrite: bool) -> tuple[str, str]:
    prefix = (case_dir / case.stem).resolve()
    nc = prefix.with_suffix(".nc")
    if result_is_complete(nc, prefix.with_suffix(".mbd")) and not overwrite:
        return "reused", "existing NetCDF retained"
    if overwrite:
        for suffix in (".nc", ".out", ".log", ".stdout", ".mov"):
            path = prefix.with_suffix(suffix)
            if path.is_file():
                path.unlink()
    result = subprocess.run(
        [str(MODEL["mbdyn_executable"]), "-s", "-f", str(prefix.with_suffix(".mbd")),
         "-o", str(prefix)],
        cwd=BASELINE_DIR,
        text=True,
        capture_output=True,
    )
    prefix.with_suffix(".stdout").write_text(result.stdout + result.stderr)
    if result.returncode != 0 or not nc.is_file():
        return "failed", f"returncode={result.returncode}; inspect .stdout"
    return "complete", "MBDyn completed"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--campaign", choices=tuple(CONFIG["campaigns"]), default="primary"
    )
    parser.add_argument(
        "--output", type=Path,
        help=f"default: {DEFAULT_OUTPUT_ROOT}/<campaign>",
    )
    parser.add_argument(
        "--execute", action="store_true",
        help="explicit opt-in; without it only .mbd inputs and manifests are prepared",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--analyse", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs deve essere almeno 1")
    output = (args.output or DEFAULT_OUTPUT_ROOT / args.campaign).expanduser().resolve()
    cases, rows, manifest = prepare(args.campaign, output, overwrite=args.overwrite)
    print(f"[prepared] {len(cases)} input in {output / 'cases'}")
    print(f"[manifest] {manifest}")
    if not args.execute:
        print("[no-run] nessuna simulazione avviata; usare --execute")
        return
    case_dir = output / "cases"
    fingerprints = source_fingerprints()
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_one, case, case_dir, args.overwrite): index
            for index, case in enumerate(cases)
        }
        for done, future in enumerate(as_completed(futures), start=1):
            index = futures[future]
            status, message = future.result()
            rows[index]["status"] = status
            rows[index]["message"] = message
            write_manifest(rows, manifest, fingerprints)
            print(f"[{done}/{len(cases)}] {rows[index]['stem']}: {status}", flush=True)
    failures = [row for row in rows if row["status"] == "failed"]
    print(f"[done] failed={len(failures)}")
    if failures:
        raise SystemExit(2)
    if args.analyse:
        command = ["python3", str(ROOT / "analyse_sweep.py"), str(output)]
        for reference in CONFIG["campaigns"][args.campaign].get(
            "reference_campaigns", []
        ):
            command.extend(["--reference-directory", str(output.parent / reference)])
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
