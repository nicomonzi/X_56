#!/usr/bin/env python3
"""Project a prestressed modal FEM snapshot onto the fixed MBDyn Phi0 basis."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


@dataclass
class FemData:
    path: Path
    node_ids: np.ndarray
    modes: np.ndarray  # (number of modes, number of nodes, 6)
    modal_mass: np.ndarray
    modal_stiffness: np.ndarray
    lumped_mass: np.ndarray


def _numbers(line: str) -> list[float]:
    return [float(item.replace("D", "E")) for item in line.split()]


def read_modal_fem(path: Path) -> FemData:
    path = path.expanduser().resolve()
    node_count = mode_count = None
    group = None
    node_ids: list[int] = []
    modes = None
    mode_index = -1
    mode_row = 0
    modal_mass: list[list[float]] = []
    modal_stiffness: list[list[float]] = []
    lumped_mass: list[list[float]] = []
    marker = re.compile(r"^\*\* RECORD GROUP (\d+),")
    mode_marker = re.compile(r"^\*\*\s+NORMAL MODE SHAPE #\s*(\d+)")
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            match = marker.match(line)
            if match:
                group = int(match.group(1))
                continue
            if line.startswith("**"):
                if group == 8:
                    match = mode_marker.match(line)
                    if match:
                        mode_index = int(match.group(1)) - 1
                        mode_row = 0
                continue
            if not line.strip():
                continue
            if group == 1 and node_count is None:
                fields = line.split()
                if len(fields) < 3:
                    raise ValueError(f"{path}:{line_number}: header FEM non valido")
                node_count, mode_count = int(fields[1]), int(fields[2])
                modes = np.empty((mode_count, node_count, 6), dtype=float)
            elif group == 2:
                node_ids.extend(int(item) for item in line.split())
            elif group == 8:
                if modes is None or mode_index < 0:
                    raise ValueError(f"{path}:{line_number}: mode shape prima dell'header")
                values = _numbers(line)
                if len(values) != 6 or mode_row >= modes.shape[1]:
                    raise ValueError(f"{path}:{line_number}: riga mode shape non valida")
                modes[mode_index, mode_row, :] = values
                mode_row += 1
            elif group == 9:
                modal_mass.append(_numbers(line))
            elif group == 10:
                modal_stiffness.append(_numbers(line))
            elif group == 11:
                lumped_mass.append(_numbers(line))
    if node_count is None or mode_count is None or modes is None:
        raise ValueError(f"{path}: header incompleto")
    if len(node_ids) != node_count:
        raise ValueError(f"{path}: attesi {node_count} nodi, trovati {len(node_ids)}")
    arrays = {
        "modal mass": np.asarray(modal_mass),
        "modal stiffness": np.asarray(modal_stiffness),
        "lumped mass": np.asarray(lumped_mass),
    }
    if arrays["modal mass"].shape != (mode_count, mode_count):
        raise ValueError(f"{path}: modal mass {arrays['modal mass'].shape}")
    if arrays["modal stiffness"].shape != (mode_count, mode_count):
        raise ValueError(f"{path}: modal stiffness {arrays['modal stiffness'].shape}")
    if arrays["lumped mass"].shape != (node_count, 6):
        raise ValueError(f"{path}: lumped mass {arrays['lumped mass'].shape}")
    return FemData(
        path, np.asarray(node_ids, dtype=int), modes,
        arrays["modal mass"], arrays["modal stiffness"], arrays["lumped mass"],
    )


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project(
    baseline: FemData,
    snapshot: FemData,
    selected_modes: list[int],
    excluded_grids: set[int],
) -> tuple[np.ndarray, dict]:
    baseline_lookup = {int(grid): i for i, grid in enumerate(baseline.node_ids)}
    snapshot_lookup = {int(grid): i for i, grid in enumerate(snapshot.node_ids)}
    common = [
        int(grid) for grid in baseline.node_ids
        if int(grid) in snapshot_lookup and int(grid) not in excluded_grids
    ]
    if not common:
        raise ValueError("nessun GRID comune tra baseline e snapshot")
    bi = np.asarray([baseline_lookup[grid] for grid in common])
    si = np.asarray([snapshot_lookup[grid] for grid in common])
    phi0 = baseline.modes[np.asarray(selected_modes) - 1][:, bi, :].transpose(1, 2, 0).reshape(-1, len(selected_modes))
    psi = snapshot.modes[:, si, :].transpose(1, 2, 0).reshape(-1, snapshot.modes.shape[0])
    weights = np.abs(baseline.lumped_mass[bi, :]).reshape(-1)
    threshold = max(float(weights.max()) * 1e-14, np.finfo(float).tiny)
    active = weights > threshold
    if np.count_nonzero(active) <= snapshot.modes.shape[0]:
        raise ValueError("DOF pesati insufficienti per il cambio di base")
    root_weight = np.sqrt(weights[active] / weights[active].max())
    a = psi[active, :] * root_weight[:, None]
    b = phi0[active, :] * root_weight[:, None]
    transform, _, rank, singular_values = np.linalg.lstsq(a, b, rcond=None)
    reconstruction = a @ transform
    residual_by_mode = (
        np.linalg.norm(reconstruction - b, axis=0)
        / np.maximum(np.linalg.norm(b, axis=0), np.finfo(float).tiny)
    )
    kh = transform.T @ snapshot.modal_stiffness @ transform
    mh = transform.T @ snapshot.modal_mass @ transform
    report = {
        "common_grid_count": len(common),
        "excluded_grids": sorted(excluded_grids),
        "weighted_dof_count": int(np.count_nonzero(active)),
        "snapshot_mode_count": int(snapshot.modes.shape[0]),
        "least_squares_rank": int(rank),
        "condition_number": float(singular_values[0] / singular_values[-1]),
        "weighted_relative_residual_by_phi0_mode": residual_by_mode.tolist(),
        "maximum_weighted_relative_residual": float(residual_by_mode.max()),
        "projected_modal_mass_identity_error": float(
            np.linalg.norm(mh - np.eye(len(selected_modes))) / np.sqrt(len(selected_modes))
        ),
        "transform": transform.tolist(),
        "projected_modal_mass": mh.tolist(),
    }
    return 0.5 * (kh + kh.T), report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-fem", type=Path, required=True)
    parser.add_argument("--load-factor", type=float, required=True, help="coordinata nominale del fit")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=ROOT / "config.json")
    parser.add_argument("--exclude-grid", type=int, action="append", default=[10062])
    parser.add_argument("--maximum-residual", type=float, default=0.05)
    parser.add_argument("--allow-large-residual", action="store_true")
    args = parser.parse_args()
    config_path = args.config.expanduser().resolve()
    config = json.loads(config_path.read_text())
    baseline_path = Path(config["model"]["baseline_fem"])
    if not baseline_path.is_absolute():
        baseline_path = (config_path.parent / baseline_path).resolve()
    baseline = read_modal_fem(baseline_path)
    snapshot = read_modal_fem(args.case_fem)
    kh, report = project(
        baseline, snapshot, list(config["model"]["selected_modes"]), set(args.exclude_grid),
    )
    if report["maximum_weighted_relative_residual"] > args.maximum_residual and not args.allow_large_residual:
        parser.error(
            f"residuo di sottospazio {report['maximum_weighted_relative_residual']:.3e} "
            f"> {args.maximum_residual:.3e}; aumentare i modi o verificare la base"
        )
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(output / "kh_fixed_basis.csv", kh, delimiter=",", fmt="%.16e")
    metadata = {
        "load_factor": args.load_factor,
        "basis_id": config["model"]["basis_id"],
        "modes": config["model"]["selected_modes"],
        "fem_sha256": config["model"]["baseline_fem_sha256"],
        "units": "MBDyn FEM modal normalization; K_h in s^-2",
        "projection": "Phi0 ~= Psi_n*T weighted subspace reconstruction; K_h= T.T@KHH_n@T",
        "case_fem": str(snapshot.path),
        "case_fem_sha256": file_hash(snapshot.path),
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (output / "projection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"output": str(output), **{k: report[k] for k in (
        "maximum_weighted_relative_residual", "projected_modal_mass_identity_error", "condition_number"
    )}}, indent=2))


if __name__ == "__main__":
    main()
