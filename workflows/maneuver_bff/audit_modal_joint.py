#!/usr/bin/env python3
"""Audit the exact MBDyn modal-joint path used by the X-56 campaign."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

import numpy as np

from campaign import CONFIG, MODEL, ROOT


def record_groups(text: str) -> list[int]:
    return [int(value) for value in re.findall(r"\*\*\s*RECORD GROUP\s+(\d+)", text)]


def fem_header(text: str) -> dict:
    match = re.search(
        r"\*\*\s*RECORD GROUP\s+1[^\n]*\n(?:\*\*[^\n]*\n)*\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)",
        text,
    )
    if not match:
        raise RuntimeError("header RECORD GROUP 1 non riconosciuto")
    return {
        "revision": match.group(1),
        "fem_nodes": int(match.group(2)),
        "normal_modes": int(match.group(3)),
        "attachment_modes": int(match.group(4)),
        "constraint_modes": int(match.group(5)),
        "rejected_modes": int(match.group(6)),
    }


def matrix_group(text: str, group: int, size: int) -> np.ndarray:
    match = re.search(
        rf"\*\*\s*RECORD GROUP\s+{group}[^\n]*\n(.*?)(?=\n\*\*\s*RECORD GROUP\s+\d+|\Z)",
        text,
        flags=re.S,
    )
    if not match:
        raise RuntimeError(f"RECORD GROUP {group} assente")
    values = [float(value.replace("D", "E")) for value in re.findall(
        r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?", match.group(1)
    )]
    expected = size * size
    if len(values) != expected:
        raise RuntimeError(
            f"RECORD GROUP {group}: attesi {expected} valori, trovati {len(values)}"
        )
    return np.asarray(values).reshape((size, size), order="F")


def selected_modes(modal_joint_text: str) -> list[int]:
    match = re.search(r"\b\d+\s*,\s*\n\s*list\s*,\s*([^,;]+(?:,[^;]+)?)", modal_joint_text)
    # The expression above is intentionally followed by a stricter match for
    # this model's six-mode list, so a later unrelated list cannot be accepted.
    strict = re.search(r"\b6\s*,\s*\n\s*list\s*,\s*([0-9,\s]+),\s*\n", modal_joint_text)
    if strict:
        return [int(value) for value in re.findall(r"\d+", strict.group(1))]
    if match:
        return [int(value) for value in re.findall(r"\d+", match.group(1))]
    raise RuntimeError("lista dei modi selezionati non riconosciuta")


def git_revision(path: Path) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        text=True,
        capture_output=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit() -> dict:
    fem_path = Path(MODEL["fem_file"])
    joint_path = Path(MODEL["modal_joint_file"])
    main_path = Path(MODEL["main_mbd_file"])
    femgen_path = Path(MODEL["femgen_source"])
    modal_source = Path(MODEL["mbdyn_modal_source"])
    modal_ad_source = Path(MODEL["mbdyn_modal_ad_source"])
    fem = fem_path.read_text(errors="replace")
    joint = joint_path.read_text()
    main = main_path.read_text()
    femgen = femgen_path.read_text(errors="replace")
    source = modal_source.read_text(errors="replace")
    source_ad = modal_ad_source.read_text(errors="replace")
    header = fem_header(fem)
    groups = record_groups(fem)
    modes = selected_modes(joint)
    mass = matrix_group(fem, 9, header["normal_modes"])
    stiffness = matrix_group(fem, 10, header["normal_modes"])
    frequencies = []
    for mode in modes:
        index = mode - 1
        frequencies.append({
            "fem_mode": mode,
            "modal_mass": float(mass[index, index]),
            "modal_stiffness": float(stiffness[index, index]),
            "frequency_hz": float(
                math.sqrt(stiffness[index, index] / mass[index, index]) / (2.0 * math.pi)
            ),
        })
    mass_offdiag = mass - np.diag(np.diag(mass))
    stiffness_offdiag = stiffness - np.diag(np.diag(stiffness))
    has_group19 = 19 in groups
    automatic_differentiation = bool(re.search(
        r"(?im)^\s*use\s+automatic\s+differentiation\s*;", main
    ))
    femgen_writes_group19 = "RECORD GROUP 19" in femgen
    source_reads_group19 = "RECORD GROUP 19" in source
    source_applies_dynamic_kgeo = all(token in source_ad for token in (
        "rgModalStressStiff", "oStressStiffIndexW", "oStressStiffIndexWP",
        "oStressStiffIndexVP", "oStressStiffIndexF", "oStressStiffIndexM",
    ))
    current_status = (
        "stress_stiffening_active"
        if has_group19 and automatic_differentiation
        else "linear_reduced_stiffness_only"
    )
    return {
        "audited_files": {
            "fem": {"path": str(fem_path), "sha256": sha256(fem_path)},
            "modal_joint": {"path": str(joint_path), "sha256": sha256(joint_path)},
            "main_mbd": {"path": str(main_path), "sha256": sha256(main_path)},
            "femgen_source": str(femgen_path),
            "modal_source": str(modal_source),
            "modal_ad_source": str(modal_ad_source),
        },
        "mbdyn_source_git_revision": git_revision(modal_source.parents[2]),
        "fem_header": header,
        "record_groups_present": groups,
        "record_group_11_lumped_mass_present": 11 in groups,
        "record_group_19_stress_stiffening_present": has_group19,
        "automatic_differentiation_enabled_in_current_model": automatic_differentiation,
        "stock_femgen_writes_record_group_19": femgen_writes_group19,
        "mbdyn_source_reads_record_group_19": source_reads_group19,
        "mbdyn_modal_ad_applies_dynamic_kgeo": source_applies_dynamic_kgeo,
        "selected_modes": frequencies,
        "mass_max_abs_offdiagonal": float(np.max(np.abs(mass_offdiag))),
        "stiffness_max_abs_offdiagonal": float(np.max(np.abs(stiffness_offdiag))),
        "current_structural_stiffness_status": current_status,
        "physical_conclusion": (
            "The current modal joint includes finite-body inertial coupling from lumped mass data, "
            "but its reduced elastic stiffness Kqq is constant. MBDyn can update Kgeo during the "
            "maneuver only when unit-load stress-stiffening matrices are supplied in RECORD GROUP 19 "
            "and automatic differentiation is enabled. The installed femgen does not generate that group."
        ),
        "campaign_decision": (
            "Run the linear-ROM refined sweep and the explicit mode-7 stiffness sensitivity first. "
            "Do not label the sensitivity force as prestress. Generate a physical prestress model only "
            "if the measured sensitivity exceeds the identification resolution."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "audit" / "modal_joint_audit.json")
    args = parser.parse_args()
    result = audit()
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2))
    print(json.dumps({
        "status": result["current_structural_stiffness_status"],
        "record_groups": result["record_groups_present"],
        "record_group_19": result["record_group_19_stress_stiffening_present"],
        "automatic_differentiation": result["automatic_differentiation_enabled_in_current_model"],
        "femgen_group_19": result["stock_femgen_writes_record_group_19"],
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
