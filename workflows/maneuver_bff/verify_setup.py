#!/usr/bin/env python3
"""Static preflight for configuration, prepared inputs, and paired coverage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from audit_modal_joint import audit
from campaign import CONFIG, ROOT, build_cases, source_fingerprints

META_RE = re.compile(r"(?m)^# STIFFNESS_STUDY_METADATA (\{.*\})$")


def verify_campaign(name: str, directory: Path) -> dict:
    expected = build_cases(name)
    manifest = directory / "manifest.csv"
    if not manifest.is_file():
        raise RuntimeError(f"manifest assente: {manifest}")
    with manifest.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    manifest_json = json.loads(manifest.with_suffix(".json").read_text())
    expected_stems = {case.stem for case in expected}
    manifest_stems = {row["stem"] for row in rows}
    inputs = {path.stem: path for path in (directory / "cases").glob("*.mbd")}
    errors: list[str] = []
    integration_steps = 0
    if len(rows) != len(expected) or manifest_stems != expected_stems:
        errors.append("manifest and configured matrix differ")
    if set(inputs) != expected_stems:
        errors.append("prepared .mbd inputs and configured matrix differ")
    if manifest_json.get("source_sha256") != source_fingerprints():
        errors.append("source fingerprints differ from the manifest")
    rows_by_stem = {row["stem"]: row for row in rows}
    pair_counts: Counter[tuple] = Counter()
    for stem, path in inputs.items():
        text = path.read_text()
        if rows_by_stem.get(stem, {}).get("input_sha256") != hashlib.sha256(
            text.encode()
        ).hexdigest():
            errors.append(f"input hash mismatch: {stem}")
        match = META_RE.search(text)
        if not match:
            errors.append(f"metadata missing: {stem}")
            continue
        meta = json.loads(match.group(1))
        key = (
            meta["velocity_mps"], meta["nominal_load_factor"], meta["time_step_s"],
            meta["mode7_frequency_scale"],
        )
        pair_counts[key] += 1
        target_dt = float(meta["time_step_s"])
        rendered_dt = re.search(r"(?m)^set: const real TIME_STEP = ([^;]+);", text)
        if not rendered_dt or abs(float(rendered_dt.group(1)) - target_dt) > 1e-14:
            errors.append(f"TIME_STEP mismatch: {stem}")
        modified = abs(float(meta["mode7_frequency_scale"]) - 1.0) > 1e-12
        has_force = "force: STIFFNESS_SCREEN_FORCE, modal" in text
        force_count = re.search(r"(?m)^\s*forces:\s*(\d+)\s*;", text)
        if has_force != modified or not force_count:
            errors.append(f"stiffness force mismatch: {stem}")
        elif int(force_count.group(1)) != (3 if modified else 2):
            errors.append(f"force count mismatch: {stem}")
        if text.count("MANEUVER_CONTROL: all aerodynamic surfaces held at release") != 1:
            errors.append(f"surface hold invariant missing: {stem}")
        final_time = re.search(r"(?m)^set: const real FINAL_TIME = ([^;]+);", text)
        release = re.search(r"(?m)^set: const real SAS_OFF_START = ([^;]+);", text)
        duration = re.search(r"(?m)^set: const real SAS_OFF_DURATION = ([^;]+);", text)
        expected_end = (
            float(release.group(1)) + float(duration.group(1))
            + int(CONFIG["execution"]["post_sas_off_margin_steps"]) * target_dt
            if release and duration else None
        )
        if (
            not final_time or expected_end is None
            or abs(float(final_time.group(1)) - expected_end) > 1e-10
        ):
            errors.append(f"optimized FINAL_TIME mismatch: {stem}")
        elif target_dt > 0.0:
            integration_steps += round(float(final_time.group(1)) / target_dt)
        for marker in (
            "default output: none;", "output meter: closest next,",
            "output: joint, 5;", "output: structural, 990000;",
        ):
            if text.count(marker) != 1:
                errors.append(f"optimized output marker mismatch ({marker}): {stem}")
    bad_pairs = [key for key, count in pair_counts.items() if count != 2]
    if bad_pairs:
        errors.append(f"unpaired conditions: {bad_pairs}")
    return {
        "campaign": name,
        "expected_inputs": len(expected),
        "prepared_inputs": len(inputs),
        "paired_conditions": len(pair_counts),
        "integration_steps": integration_steps,
        "errors": errors,
        "ok": not errors,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs", type=Path, default=Path(CONFIG["execution"]["output_root"]),
    )
    parser.add_argument("--output", type=Path, default=ROOT / "audit" / "preflight.json")
    args = parser.parse_args()
    reports = [
        verify_campaign(name, args.runs / name) for name in CONFIG["campaigns"]
    ]
    modal = audit()
    optimized_trajectories = sum(item["expected_inputs"] for item in reports)
    optimized_steps = sum(item["integration_steps"] for item in reports)
    runtime_reference = CONFIG["execution"]["runtime_reference"]
    seconds_per_step = float(runtime_reference["median_wall_seconds_per_integration_step"])
    old_trajectories = int(runtime_reference["superseded_trajectory_count"])
    old_steps = int(runtime_reference["superseded_integration_steps"])
    report = {
        "campaigns": reports,
        "optimization": {
            "optimized_trajectories": optimized_trajectories,
            "superseded_trajectories": old_trajectories,
            "trajectory_reduction_percent": 100.0 * (1.0 - optimized_trajectories / old_trajectories),
            "optimized_integration_steps": optimized_steps,
            "superseded_integration_steps": old_steps,
            "integration_step_reduction_percent": 100.0 * (1.0 - optimized_steps / old_steps),
            "estimated_serial_wall_hours": optimized_steps * seconds_per_step / 3600.0,
            "runtime_reference": runtime_reference,
        },
        "source_sha256": source_fingerprints(),
        "modal_joint_status": modal["current_structural_stiffness_status"],
        "record_group_19_present": modal["record_group_19_stress_stiffening_present"],
        "all_checks_passed": all(item["ok"] for item in reports),
    }
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if not report["all_checks_passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
