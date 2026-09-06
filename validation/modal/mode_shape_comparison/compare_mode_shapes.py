#!/usr/bin/env python3
"""Compare Nastran F06 eigenvectors with the shapes used by an MBDyn modal joint."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_FEM = Path("../NASTRAN/FEMGEN40/mbdyn_modal.fem")
DEFAULT_F06 = Path("nastran/sol103_10lb_f06.f06")
SPAN_NODE_IDS = list(range(990001, 990024)) + list(range(991002, 991024))
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EDed][-+]?\d+)?")


def nastran_float(token: str) -> float:
    """Read ordinary, D-exponent, and Nastran compact-exponent numbers."""
    value = token.strip().replace("D", "E").replace("d", "e")
    if "e" not in value.lower():
        value = re.sub(r"(?<=\d)([+-]\d+)$", r"E\1", value)
    return float(value)


def numbers(line: str) -> list[float]:
    return [nastran_float(token) for token in FLOAT_RE.findall(line)]


def seek_marker(lines, marker: str) -> None:
    for line in lines:
        if marker in line:
            return
    raise ValueError(f"Missing FEM marker: {marker}")


def values_until(lines, marker: str, integer: bool = False) -> list[float] | list[int]:
    values = []
    for line in lines:
        if marker in line:
            return values
        if line.lstrip().startswith("**"):
            continue
        if integer:
            values.extend(int(token) for token in re.findall(r"\d+", line))
        else:
            values.extend(numbers(line))
    raise ValueError(f"Missing FEM marker: {marker}")


def parse_fem(path: Path, wanted_modes: set[int]) -> tuple[
    dict[int, np.ndarray], dict[int, np.ndarray], dict[int, float]
]:
    """Stream the large femgen ASCII file and retain only span-line data."""
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        lines = iter(stream)
        seek_marker(lines, "RECORD GROUP 1")
        header = None
        for line in lines:
            if line.lstrip().startswith("**"):
                continue
            fields = line.split()
            if len(fields) >= 3:
                header = fields
                break
        if header is None:
            raise ValueError("Could not read FEM header")
        node_count, mode_count = int(header[1]), int(header[2])

        seek_marker(lines, "RECORD GROUP 2")
        node_ids = values_until(lines, "RECORD GROUP 3", integer=True)
        if len(node_ids) != node_count:
            raise ValueError(f"FEM declares {node_count} nodes but lists {len(node_ids)}")
        node_index = {node_id: index for index, node_id in enumerate(node_ids)}
        missing = sorted(set(SPAN_NODE_IDS) - node_index.keys())
        if missing:
            raise ValueError(f"Span nodes missing from FEM: {missing}")

        seek_marker(lines, "RECORD GROUP 5")
        x = np.asarray(values_until(lines, "RECORD GROUP 6"), dtype=float)
        y = np.asarray(values_until(lines, "RECORD GROUP 7"), dtype=float)
        z = np.asarray(values_until(lines, "RECORD GROUP 8"), dtype=float)
        if any(len(axis) != node_count for axis in (x, y, z)):
            raise ValueError("FEM coordinate record size does not match its node count")

        coordinates = {
            node_id: np.array([x[node_index[node_id]], y[node_index[node_id]], z[node_index[node_id]]])
            for node_id in SPAN_NODE_IDS
        }
        shapes: dict[int, np.ndarray] = {}
        retained_indices = {node_index[node_id]: row for row, node_id in enumerate(SPAN_NODE_IDS)}

        modes_read = 0
        for line in lines:
            if "RECORD GROUP 9" in line:
                break
            match = re.search(r"NORMAL MODE SHAPE #\s*(\d+)", line)
            if not match:
                continue
            mode = int(match.group(1))
            modes_read += 1
            selected = np.zeros((len(SPAN_NODE_IDS), 6), dtype=float) if mode in wanted_modes else None
            for index in range(node_count):
                try:
                    row = numbers(next(lines))
                except StopIteration as exc:
                    raise ValueError(f"Unexpected EOF while reading FEM mode {mode}") from exc
                if len(row) != 6:
                    raise ValueError(f"FEM mode {mode}, node row {index + 1}: expected 6 values")
                if selected is not None and index in retained_indices:
                    selected[retained_indices[index], :] = row
            if selected is not None:
                shapes[mode] = selected

        if modes_read != mode_count:
            raise ValueError(f"FEM declares {mode_count} modes but contains {modes_read}")

        seek_marker(lines, "RECORD GROUP 10")
        stiffness_values = np.asarray(values_until(lines, "RECORD GROUP 11"), dtype=float)
        if stiffness_values.size != mode_count * mode_count:
            raise ValueError("FEM modal stiffness matrix has an unexpected size")
        stiffness = stiffness_values.reshape(mode_count, mode_count)
        frequencies = {
            mode: math.sqrt(max(stiffness[mode - 1, mode - 1], 0.0)) / (2.0 * math.pi)
            for mode in wanted_modes
            if 1 <= mode <= mode_count
        }

    return coordinates, shapes, frequencies


def parse_f06(path: Path) -> tuple[dict[int, dict[int, np.ndarray]], dict[int, float]]:
    """Read printed modal vectors and the eigenvalue table from an MSC Nastran F06."""
    shapes: dict[int, dict[int, np.ndarray]] = {}
    frequencies: dict[int, float] = {}
    current_mode: int | None = None
    in_eigenvalue_table = False

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            compact = re.sub(r"\s+", "", line).upper()
            if "REALEIGENVALUES" in compact:
                in_eigenvalue_table = True
                current_mode = None
                continue

            if in_eigenvalue_table:
                fields = line.split()
                if len(fields) >= 7 and fields[0].isdigit() and fields[1].isdigit():
                    try:
                        mode = int(fields[0])
                        frequencies[mode] = nastran_float(fields[4])
                        continue
                    except ValueError:
                        pass
                if "USER INFORMATION" in line or "EIGENVALUE =" in line:
                    in_eigenvalue_table = False

            mode_match = re.search(r"EIGENVECTORNO\.?(\d+)", compact)
            if mode_match:
                current_mode = int(mode_match.group(1))
                shapes.setdefault(current_mode, {})
                continue

            if current_mode is None:
                continue
            fields = line.split()
            if len(fields) < 8 or not fields[0].isdigit() or fields[1].upper() not in {"G", "S"}:
                continue
            node_id = int(fields[0])
            if node_id not in SPAN_NODE_IDS:
                continue
            try:
                shapes[current_mode][node_id] = np.asarray(
                    [nastran_float(value) for value in fields[2:8]], dtype=float
                )
            except ValueError:
                continue

    shapes = {
        mode: node_data
        for mode, node_data in shapes.items()
        if set(SPAN_NODE_IDS).issubset(node_data)
    }
    return shapes, frequencies


def parse_mbdyn_mov(path: Path) -> np.ndarray:
    """Extract the largest physical displacement snapshot from an MBDyn MOV file."""
    snapshots: list[dict[int, np.ndarray]] = []
    current: dict[int, np.ndarray] = {}
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 4 or not fields[0].isdigit():
                continue
            node_id = int(fields[0])
            if node_id not in SPAN_NODE_IDS:
                continue
            current[node_id] = np.asarray([float(value) for value in fields[1:4]])
            if len(current) == len(SPAN_NODE_IDS):
                snapshots.append(current)
                current = {}
    if len(snapshots) < 2:
        raise ValueError("MBDyn MOV must contain an undeformed and at least one deformed snapshot")
    reference = np.vstack([snapshots[0][node_id] for node_id in SPAN_NODE_IDS])
    displacements = [
        np.vstack([snapshot[node_id] for node_id in SPAN_NODE_IDS]) - reference
        for snapshot in snapshots[1:]
    ]
    selected = max(displacements, key=lambda value: np.linalg.norm(value))
    result = np.zeros((len(SPAN_NODE_IDS), 6), dtype=float)
    result[:, :3] = selected
    return result


def parse_modes(specification: str) -> list[int]:
    modes: set[int] = set()
    for item in specification.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            start, stop = (int(value) for value in item.split("-", maxsplit=1))
            modes.update(range(start, stop + 1))
        else:
            modes.add(int(item))
    return sorted(modes)


def spanwise_arc_length(coordinates: dict[int, np.ndarray]) -> np.ndarray:
    """Return signed elastic-axis arc length, including the swept/vertical winglet."""
    station: dict[int, float] = {990001: 0.0}
    previous = coordinates[990001]
    distance = 0.0
    for node_id in range(990002, 990024):
        distance += float(np.linalg.norm(coordinates[node_id] - previous))
        station[node_id] = -distance
        previous = coordinates[node_id]

    previous = coordinates[990001]
    distance = 0.0
    for node_id in range(991002, 991024):
        distance += float(np.linalg.norm(coordinates[node_id] - previous))
        station[node_id] = distance
        previous = coordinates[node_id]
    return np.asarray([station[node_id] for node_id in SPAN_NODE_IDS])


def normalized_and_aligned(nastran: np.ndarray, mbdyn: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nastran_xyz = nastran[:, :3].copy()
    mbdyn_xyz = mbdyn[:, :3].copy()
    nastran_scale = np.max(np.linalg.norm(nastran_xyz, axis=1))
    mbdyn_scale = np.max(np.linalg.norm(mbdyn_xyz, axis=1))
    if nastran_scale <= 0.0 or mbdyn_scale <= 0.0:
        raise ValueError("A mode shape has zero translational amplitude on all comparison nodes")
    nastran_xyz /= nastran_scale
    mbdyn_xyz /= mbdyn_scale
    if np.vdot(nastran_xyz.ravel(), mbdyn_xyz.ravel()).real < 0.0:
        mbdyn_xyz *= -1.0
    return nastran_xyz, mbdyn_xyz


def modal_assurance_criterion(left: np.ndarray, right: np.ndarray) -> float:
    a, b = left.ravel(), right.ravel()
    denominator = float(np.vdot(a, a).real * np.vdot(b, b).real)
    return float(abs(np.vdot(a, b)) ** 2 / denominator) if denominator else float("nan")


def compare(args: argparse.Namespace) -> int:
    modes = parse_modes(args.modes)
    fem_path = args.fem.resolve()
    f06_path = args.f06.resolve()
    output_dir = args.output.resolve()

    coordinates, fem_shapes, fem_frequencies = parse_fem(fem_path, set(modes))
    f06_shapes, f06_frequencies = parse_f06(f06_path)
    mov_mode = None
    if args.mbdyn_mov is not None:
        if args.mov_mode is not None:
            mov_mode = args.mov_mode
        elif len(modes) == 1:
            mov_mode = modes[0]
        else:
            raise ValueError("--mov-mode is required when --mbdyn-mov is used with multiple modes")
        if mov_mode not in modes:
            raise ValueError("--mov-mode must be included in --modes")
        fem_shapes[mov_mode] = parse_mbdyn_mov(args.mbdyn_mov.resolve())

    available = [mode for mode in modes if mode in fem_shapes and mode in f06_shapes]
    missing = [mode for mode in modes if mode not in available]
    if not available:
        raise ValueError(
            "No complete printed F06 eigenvectors were found. Run the supplied deck; "
            "the old NASTRAN40 F06 used PLOT-only output and cannot be parsed for shapes."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    span = spanwise_arc_length(coordinates)
    order = np.argsort(span)
    rows = []

    for mode in available:
        nastran_raw = np.vstack([f06_shapes[mode][node_id] for node_id in SPAN_NODE_IDS])
        nastran, mbdyn = normalized_and_aligned(nastran_raw, fem_shapes[mode])
        difference = nastran - mbdyn
        mac = modal_assurance_criterion(nastran, mbdyn)
        rms = float(np.sqrt(np.mean(difference**2)))
        max_error = float(np.max(np.linalg.norm(difference, axis=1)))
        f_nas = f06_frequencies.get(mode, float("nan"))
        f_mbd = fem_frequencies.get(mode, float("nan"))
        freq_error = 100.0 * (f_mbd - f_nas) / f_nas if f_nas else float("nan")
        rows.append([mode, f_nas, f_mbd, freq_error, mac, rms, max_error])

        fig, axes = plt.subplots(3, 1, figsize=(9.0, 10.0), sharex=True)
        labels = ("Normalized X displacement", "Normalized Y displacement", "Normalized Z displacement")
        for component, (axis, label) in enumerate(zip(axes, labels)):
            axis.plot(span[order], nastran[order, component], "o-", label="Nastran SOL 103 (F06)")
            mbdyn_label = "MBDyn modal-joint run (.mov)" if mode == mov_mode else "MBDyn modal joint (.fem)"
            axis.plot(span[order], mbdyn[order, component], "s--", label=mbdyn_label)
            axis.plot(span[order], difference[order, component], ":", color="0.35", label="Difference")
            axis.set_ylabel(label)
            axis.grid(True, alpha=0.3)
        axes[0].legend(loc="best", fontsize=8)
        axes[-1].set_xlabel("Signed spanwise arc length [in]")
        fig.suptitle(
            f"Mode {mode}: Nastran {f_nas:.4f} Hz, MBDyn {f_mbd:.4f} Hz, MAC {mac:.6f}"
        )
        fig.tight_layout()
        fig.savefig(output_dir / f"mode_{mode:02d}_spanwise_comparison.png", dpi=180)
        plt.close(fig)

    with (output_dir / "modal_comparison_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "mode",
                "nastran_frequency_hz",
                "mbdyn_fem_frequency_hz",
                "frequency_difference_percent",
                "translational_MAC",
                "normalized_RMS_difference",
                "maximum_nodal_vector_difference",
            ]
        )
        writer.writerows(rows)

    print(f"Compared {len(available)} modes. Results: {output_dir}")
    if missing:
        print(f"Skipped modes without complete data: {', '.join(map(str, missing))}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--f06", type=Path, default=base / DEFAULT_F06, help="Nastran F06 with printed eigenvectors")
    parser.add_argument("--fem", type=Path, default=(base / DEFAULT_FEM), help="femgen file used by MBDyn")
    parser.add_argument("--modes", default="7-18", help="Mode list/ranges, e.g. 7-18 or 7,9,12")
    parser.add_argument("--mbdyn-mov", type=Path, help="Optional MBDyn MOV file for a direct run comparison")
    parser.add_argument("--mov-mode", type=int, help="Mode represented by --mbdyn-mov")
    parser.add_argument("--output", type=Path, default=base / "plots", help="Output directory")
    return parser


if __name__ == "__main__":
    try:
        raise SystemExit(compare(build_parser().parse_args()))
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
