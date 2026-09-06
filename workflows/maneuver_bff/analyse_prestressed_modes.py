#!/usr/bin/env python3
"""Match prestressed SOL103 modes to the validated baseline on 45 interface grids."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[1]
BASELINE_F06 = REPO_ROOT / "validation/modal/sol103_60_modes/MAIN/sol103_60_modes.f06"
BASELINE_FEM = REPO_ROOT / "assets/fem/mbdyn_modal_60.fem"
DEFAULT_RESULTS = Path("/mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/load_recovery/prestress_loads")
SPAN_NODES = list(range(990001, 990024)) + list(range(991002, 991024))


def nastran_float(value: str) -> float:
    value = value.replace("D", "E").replace("d", "E")
    if "E" not in value.upper():
        value = re.sub(r"(?<=\d)([+-]\d+)$", r"E\1", value)
    return float(value)


def parse_f06(
    path: Path, require_shapes: bool = True
) -> tuple[dict[int, float], dict[int, np.ndarray], dict[int, float]]:
    frequencies: dict[int, float] = {}
    eigenvalues: dict[int, float] = {}
    values: dict[int, dict[int, np.ndarray]] = {}
    current_mode: int | None = None
    in_eigenvalues = False
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "FATAL MESSAGE" in line:
                raise RuntimeError(f"{path}: {line.strip()}")
            compact = re.sub(r"\s+", "", line).upper()
            if "REALEIGENVALUES" in compact:
                in_eigenvalues = True
                current_mode = None
                continue
            if in_eigenvalues:
                fields = line.split()
                if len(fields) >= 7 and fields[0].isdigit() and fields[1].isdigit():
                    try:
                        eigenvalues[int(fields[0])] = nastran_float(fields[2])
                        frequencies[int(fields[0])] = nastran_float(fields[4])
                        continue
                    except ValueError:
                        pass
                if "EIGENVALUE=" in compact or "USERINFORMATIONMESSAGE" in compact:
                    in_eigenvalues = False
            match = re.search(r"EIGENVECTORNO\.?([0-9]+)", compact)
            if match:
                current_mode = int(match.group(1))
                values.setdefault(current_mode, {})
                continue
            if current_mode is None:
                continue
            fields = line.split()
            if len(fields) >= 8 and fields[0].isdigit() and fields[1].upper() == "G":
                node = int(fields[0])
                if node in SPAN_NODES:
                    values[current_mode][node] = np.asarray(
                        [nastran_float(item) for item in fields[2:8]], dtype=float
                    )
    shapes = {
        mode: np.vstack([nodes[node] for node in SPAN_NODES])
        for mode, nodes in values.items()
        if all(node in nodes for node in SPAN_NODES)
    }
    if not frequencies or (require_shapes and not shapes):
        raise RuntimeError(f"{path}: eigenvalues or 45-grid PRINT eigenvectors missing")
    return frequencies, shapes, eigenvalues


def parse_baseline_fem(path: Path, modes: set[int]) -> dict[int, np.ndarray]:
    node_labels: list[int] = []
    shapes: dict[int, dict[int, np.ndarray]] = {}
    group = 0
    mode: int | None = None
    node_index = 0
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            marker = re.search(r"RECORD GROUP\s+(\d+)", line)
            if marker:
                group = int(marker.group(1))
                if group > 8:
                    break
                continue
            if group == 2 and not line.lstrip().startswith("*"):
                node_labels.extend(int(item) for item in line.split())
                continue
            if group != 8:
                continue
            match = re.search(r"NORMAL MODE SHAPE\s*#\s*(\d+)", line)
            if match:
                mode = int(match.group(1))
                node_index = 0
                if mode in modes:
                    shapes[mode] = {}
                continue
            if mode is None or line.lstrip().startswith("*") or not line.split():
                continue
            fields = line.split()
            if len(fields) != 6:
                continue
            if node_index >= len(node_labels):
                raise RuntimeError("FEM mode-shape block longer than node list")
            node = node_labels[node_index]
            if mode in modes and node in SPAN_NODES:
                shapes[mode][node] = np.asarray([float(item) for item in fields])
            node_index += 1
    result = {
        mode: np.vstack([values[node] for node in SPAN_NODES])
        for mode, values in shapes.items()
        if all(node in values for node in SPAN_NODES)
    }
    if set(result) != modes:
        raise RuntimeError(f"baseline FEM incomplete for modes {sorted(modes - set(result))}")
    return result


def mac(left: np.ndarray, right: np.ndarray) -> float:
    a, b = left.ravel(), right.ravel()
    denominator = float(np.vdot(a, a).real * np.vdot(b, b).real)
    return float(abs(np.vdot(a, b)) ** 2 / denominator) if denominator else 0.0


def match_case(
    baseline_frequencies: dict[int, float],
    baseline_shapes: dict[int, np.ndarray],
    f06: Path,
) -> tuple[list[dict], list[dict]]:
    frequencies, shapes, eigenvalues = parse_f06(f06)
    reference_modes = sorted(baseline_shapes)
    # The prestressed screen retains the CG support in the modal subcase, so
    # elastic roots start at 1 rather than after six free-free rigid modes.
    candidates = sorted(
        mode for mode in shapes
        if eigenvalues.get(mode, 0.0) > 0.0
        and frequencies.get(mode, 0.0) > 0.1
        and mode <= 30
    )
    mac_matrix = np.asarray([
        [mac(baseline_shapes[left], shapes[right]) for right in candidates]
        for left in reference_modes
    ])
    gaps = np.asarray([
        [abs(frequencies[right] - baseline_frequencies[left]) / baseline_frequencies[left]
         for right in candidates]
        for left in reference_modes
    ])
    cost = 1.0 - mac_matrix + 0.25 * gaps
    cost[gaps > 0.50] = 100.0 + gaps[gaps > 0.50]
    row_indices, column_indices = linear_sum_assignment(cost)
    rows = []
    for i, j in sorted(zip(row_indices, column_indices)):
        left, right = reference_modes[i], candidates[j]
        base = baseline_frequencies[left]
        current = frequencies[right]
        rows.append({
            "baseline_mode": left,
            "prestressed_mode": right,
            "baseline_frequency_hz": base,
            "prestressed_frequency_hz": current,
            "shift_from_unloaded_percent": 100.0 * (current - base) / base,
            "MAC": float(mac_matrix[i, j]),
            "reliable_MAC_ge_0p90": bool(mac_matrix[i, j] >= 0.90),
        })
    negative = [
        {"mode": mode, "eigenvalue": value, "frequency_column_hz": frequencies.get(mode)}
        for mode, value in sorted(eigenvalues.items())
        # Ignore only the numerical near-zero band below 0.1 Hz.  In the
        # CG-supported problem a more negative root is an elastic instability.
        if value < -(2.0 * np.pi * 0.1) ** 2
    ]
    return rows, negative


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--baseline-f06", type=Path, default=BASELINE_F06)
    parser.add_argument("--baseline-fem", type=Path, default=BASELINE_FEM)
    parser.add_argument("--threshold-percent", type=float, default=0.4)
    args = parser.parse_args()
    case_paths = {
        variant: {
            tag: args.results / f"{tag}/prestressed_modes_{variant}.f06"
            for tag in ("n1p0", "n1p6")
        }
        for variant in ("nodlm", "with_dlm")
    }
    supported_unloaded_path = args.results / "n1p0/supported_unloaded.f06"
    missing = [
        str(path) for paths in case_paths.values() for path in paths.values()
        if not path.is_file()
    ]
    if not supported_unloaded_path.is_file():
        missing.append(str(supported_unloaded_path))
    if missing:
        raise SystemExit("Run the four supplied Nastran decks first; missing: " + ", ".join(missing))
    baseline_frequencies, _, _ = parse_f06(args.baseline_f06, require_shapes=False)
    modes = set(range(7, 13))
    baseline_shapes = parse_baseline_fem(args.baseline_fem, modes)
    supported_rows, supported_negative = match_case(
        baseline_frequencies, baseline_shapes, supported_unloaded_path
    )
    supported_by_mode = {row["baseline_mode"]: row for row in supported_rows}
    combined = []
    mode7_by_variant = {}
    negative_eigenvalues = {}
    for variant, paths in case_paths.items():
        parsed = {
            tag: match_case(baseline_frequencies, baseline_shapes, path)
            for tag, path in paths.items()
        }
        matches = {tag: result[0] for tag, result in parsed.items()}
        negative_eigenvalues[variant] = {
            tag: result[1] for tag, result in parsed.items()
        }
        by_case = {
            tag: {row["baseline_mode"]: row for row in rows}
            for tag, rows in matches.items()
        }
        for mode in sorted(modes):
            low, high = by_case["n1p0"][mode], by_case["n1p6"][mode]
            unloaded = supported_by_mode[mode]
            incremental = 100.0 * (
                high["prestressed_frequency_hz"] - low["prestressed_frequency_hz"]
            ) / low["prestressed_frequency_hz"]
            from_unloaded_low = 100.0 * (
                low["prestressed_frequency_hz"] - unloaded["prestressed_frequency_hz"]
            ) / unloaded["prestressed_frequency_hz"]
            from_unloaded_high = 100.0 * (
                high["prestressed_frequency_hz"] - unloaded["prestressed_frequency_hz"]
            ) / unloaded["prestressed_frequency_hz"]
            row = {
                "variant": variant,
                "baseline_mode": mode,
                "n1p0_mode": low["prestressed_mode"],
                "n1p6_mode": high["prestressed_mode"],
                "f_unloaded_hz": low["baseline_frequency_hz"],
                "f_supported_unloaded_hz": unloaded["prestressed_frequency_hz"],
                "f_n1p0_hz": low["prestressed_frequency_hz"],
                "f_n1p6_hz": high["prestressed_frequency_hz"],
                "unloaded_to_n1p0_shift_percent": from_unloaded_low,
                "unloaded_to_n1p6_shift_percent": from_unloaded_high,
                "n1p0_to_n1p6_shift_percent": incremental,
                "MAC_supported_unloaded": unloaded["MAC"],
                "MAC_n1p0": low["MAC"],
                "MAC_n1p6": high["MAC"],
                "reliable": (
                    unloaded["reliable_MAC_ge_0p90"]
                    and low["reliable_MAC_ge_0p90"]
                    and high["reliable_MAC_ge_0p90"]
                ),
            }
            combined.append(row)
            if mode == 7:
                mode7_by_variant[variant] = row
    shifts = [row["unloaded_to_n1p6_shift_percent"] for row in mode7_by_variant.values()]
    envelopes = [
        max(
            abs(row["unloaded_to_n1p0_shift_percent"]),
            abs(row["unloaded_to_n1p6_shift_percent"]),
        )
        for row in mode7_by_variant.values()
    ]
    reliable = all(row["reliable"] for row in mode7_by_variant.values())
    resolved = reliable and all(value >= args.threshold_percent for value in envelopes)
    same_sign = shifts[0] * shifts[1] > 0.0
    negative_eigenvalues["supported_unloaded"] = supported_negative
    has_negative = any(
        roots
        for key, cases in negative_eigenvalues.items()
        for roots in (cases.values() if isinstance(cases, dict) else [cases])
    )
    if has_negative:
        decision = "negative_elastic_eigenvalues_prestressed_state_unstable_or_invalid"
    elif resolved and same_sign:
        decision = "physical_prestress_resolved_robust_to_dlm_distribution_build_rom"
    elif reliable and all(value < args.threshold_percent for value in envelopes):
        decision = "physical_prestress_below_resolution_skip_rom"
    else:
        decision = "physical_prestress_sensitive_to_dlm_distribution_review_before_rom"
    args.results.mkdir(parents=True, exist_ok=True)
    with (args.results / "prestressed_mode_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=combined[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(combined)
    summary = {
        "comparison_basis": (
            "frequency shifts use identical CG-supported boundary conditions; "
            "mode identity uses MAC on all 6 DOF of 45 interface grids"
        ),
        "threshold_percent": args.threshold_percent,
        "mode7_incremental_shift_percent": {
            variant: row["n1p0_to_n1p6_shift_percent"]
            for variant, row in mode7_by_variant.items()
        },
        "mode7_shift_from_supported_unloaded_percent": {
            variant: {
                "n1p0": row["unloaded_to_n1p0_shift_percent"],
                "n1p6": row["unloaded_to_n1p6_shift_percent"],
            }
            for variant, row in mode7_by_variant.items()
        },
        "negative_eigenvalues": negative_eigenvalues,
        "decision": decision,
        "modes": combined,
    }
    (args.results / "prestressed_mode_comparison.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
