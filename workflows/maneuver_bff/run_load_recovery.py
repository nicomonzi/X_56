#!/usr/bin/env python3
"""Prepare and optionally run the two shadow cases used to recover preload loads."""

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

from campaign import BASELINE_DIR, CONFIG, MODEL, StudyCase, render_case, source_fingerprints


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path(CONFIG["execution"]["output_root"]) / "load_recovery"
INTERFACE_NODES = [
    *range(990001, 990024),
    *range(991002, 991024),
]
AERODYNAMIC_ELEMENTS = list(range(1, 59))


def _constant(text: str, name: str) -> float:
    match = re.search(
        rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*([-+0-9.eE]+)\s*;",
        text,
    )
    if not match:
        raise RuntimeError(f"costante {name} non trovata")
    return float(match.group(1))


def recovery_cases() -> list[StudyCase]:
    # 66.75 m/s is the common near-onset point.  The two load states isolate
    # the increment caused by the highest accepted pull-up from steady 1 g.
    return [
        StudyCase("primary", 66.75, nominal_n, 0.01, 1.0, False)
        for nominal_n in (1.0, 1.6)
    ]


def recovery_stem(case: StudyCase) -> str:
    return case.stem.replace("primary__", "load_recovery__", 1)


def render_recovery(case: StudyCase) -> str:
    text = render_case(case)
    text = text.replace(
        '"campaign":"primary"', '"campaign":"load_recovery"', 1
    )
    text = text.replace(
        "# STIFFNESS_STUDY_METADATA ",
        "# LOAD_RECOVERY_PURPOSE body-frame averaged aerodynamic preload\n"
        "# STIFFNESS_STUDY_METADATA ",
        1,
    )

    # The modal base acceleration is needed for the inertial balance.  The
    # 45 interface nodes are needed to move each Gauss-point resultant to its
    # corresponding Nastran RBE3 reference grid.
    text, count = re.subn(
        r"(?m)^\s*output:\s*structural,\s*990000\s*;",
        "    output: structural, accelerations, 990000;",
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError("output del nodo modale 990000 non trovato")

    nodes_marker = "end: nodes;"
    if text.count(nodes_marker) != 1:
        raise RuntimeError("fine sezione nodes non univoca")
    extra_nodes = "".join(
        f"    output: structural, {label};\n"
        for label in INTERFACE_NODES
        if f"output: structural, {label};" not in text
    )
    text = text.replace(nodes_marker, extra_nodes + nodes_marker, 1)

    # Position plus force and moment at the three Gauss points is sufficient.
    # The production output meter is retained, so only the final portion of
    # the maneuver is written even though the full trajectory is integrated.
    text = text.replace(
        "default aerodynamic output: position, velocity, force, moment;",
        "default aerodynamic output: position, force, moment;",
        1,
    )
    aero_include = f'    include: "{BASELINE_DIR / "INCLUDE/aerobody.mbd"}";'
    marker = aero_include + "\n"
    if text.count(marker) != 1:
        raise RuntimeError("include aerodinamico non univoco")
    # Do not use ``range`` here: in this MBDyn build the range branch resets
    # the element flag to the legacy coefficient-only output and drops the
    # custom position/force/moment mask.  Explicit labels preserve that mask.
    aero_outputs = (
        "    output: aerodynamic body, "
        + ", ".join(str(label) for label in AERODYNAMIC_ELEMENTS)
        + ";\n"
    )
    text = text.replace(marker, marker + aero_outputs, 1)

    required = (
        "LOAD_RECOVERY_PURPOSE",
        "output: structural, accelerations, 990000",
        "output: structural, 990023",
        "output: structural, 991023",
        "output: aerodynamic body, 1, 2, 3",
        "default aerodynamic output: position, force, moment",
        "BFF_RAP_AMPLITUDE = 0.0000000000e+00",
    )
    if not all(item in text for item in required):
        raise RuntimeError("invarianti del caso load-recovery non soddisfatte")
    return text


def complete(nc: Path, mbd: Path) -> bool:
    if not nc.is_file():
        return False
    text = mbd.read_text()
    target = _constant(text, "SAS_OFF_START") + _constant(text, "SAS_OFF_DURATION")
    try:
        with Dataset(nc) as data:
            time = data["time"][:]
            return bool(
                len(time) >= 100
                and float(time[-1]) >= target - 0.011
                and "elem.aerodynamic.58.F_2" in data.variables
            )
    except Exception:
        return False


def write_manifest(rows: list[dict], output: Path) -> None:
    fields = list(rows[0])
    with (output / "manifest.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "manifest.json").write_text(json.dumps({
        "purpose": "MBDyn shadow load recovery for physical prestress",
        "source_sha256": source_fingerprints(),
        "cases": rows,
    }, indent=2))


def prepare(output: Path, overwrite: bool) -> tuple[list[StudyCase], list[dict]]:
    case_dir = output / "cases"
    case_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    cases = recovery_cases()
    for index, case in enumerate(cases, 1):
        stem = recovery_stem(case)
        mbd = case_dir / f"{stem}.mbd"
        rendered = render_recovery(case)
        nc = mbd.with_suffix(".nc")
        if mbd.is_file() and nc.is_file() and mbd.read_text() != rendered and not overwrite:
            raise RuntimeError(f"input cambiato ma risultato esistente: {nc}")
        mbd.write_text(rendered)
        is_complete = complete(nc, mbd) and not overwrite
        rows.append({
            "case_index": index,
            "velocity_mps": case.velocity_mps,
            "nominal_load_factor": case.nominal_load_factor,
            "time_step_s": case.time_step_s,
            "excited": False,
            "stem": stem,
            "input_sha256": hashlib.sha256(rendered.encode()).hexdigest(),
            "status": "existing" if is_complete else "prepared",
            "message": "complete load output retained" if is_complete else "input prepared",
        })
    write_manifest(rows, output)
    return cases, rows


def run_one(case: StudyCase, output: Path, overwrite: bool) -> tuple[str, str]:
    prefix = (output / "cases" / recovery_stem(case)).resolve()
    mbd = prefix.with_suffix(".mbd")
    nc = prefix.with_suffix(".nc")
    if complete(nc, mbd) and not overwrite:
        return "reused", "complete NetCDF retained"
    if overwrite:
        for suffix in (".nc", ".out", ".log", ".stdout", ".mov"):
            path = prefix.with_suffix(suffix)
            if path.is_file():
                path.unlink()
    result = subprocess.run(
        [str(MODEL["mbdyn_executable"]), "-s", "-f", str(mbd), "-o", str(prefix)],
        cwd=BASELINE_DIR,
        text=True,
        capture_output=True,
    )
    prefix.with_suffix(".stdout").write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        return "failed", f"returncode={result.returncode}; inspect .stdout"
    if not complete(nc, mbd):
        return "failed", "solver ended but required load variables are incomplete"
    return "complete", "MBDyn completed with aerodynamic loads; base accelerations are differentiated from velocity"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs deve essere almeno 1")
    output = args.output.expanduser().resolve()
    cases, rows = prepare(output, args.overwrite)
    print(f"[prepared] {len(cases)} shadow inputs in {output / 'cases'}")
    if not args.execute:
        print("[no-run] use --execute to run MBDyn")
        return
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_one, case, output, args.overwrite): index
            for index, case in enumerate(cases)
        }
        for done, future in enumerate(as_completed(futures), 1):
            index = futures[future]
            status, message = future.result()
            rows[index]["status"] = status
            rows[index]["message"] = message
            write_manifest(rows, output)
            print(f"[{done}/{len(cases)}] {rows[index]['stem']}: {status}", flush=True)
    failures = [row for row in rows if row["status"] == "failed"]
    print(f"[done] failed={len(failures)}")
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
