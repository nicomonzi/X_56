#!/usr/bin/env python3
"""Fit K_h,n in a fixed modal basis and generate the MBDyn modal force."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class MatrixCase:
    load_factor: float
    matrix: np.ndarray
    metadata: dict[str, Any]
    directory: Path
    symmetry_error: float


def load_config(path: Path | None = None) -> tuple[dict[str, Any], Path]:
    config_path = (path or ROOT / "config.json").expanduser().resolve()
    return json.loads(config_path.read_text()), config_path.parent


def resolve(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fem_matrix(path: Path, group: int) -> np.ndarray:
    """Read a square numeric matrix record group from an MBDyn modal FEM file."""
    lines = path.read_text().splitlines()
    marker = re.compile(rf"^\*\* RECORD GROUP {group},")
    start = next((i + 1 for i, line in enumerate(lines) if marker.match(line)), None)
    if start is None:
        raise ValueError(f"record group {group} non trovato in {path}")
    rows: list[list[float]] = []
    for line in lines[start:]:
        if line.startswith("**"):
            break
        if line.strip():
            rows.append([float(item.replace("D", "E")) for item in line.split()])
    matrix = np.asarray(rows, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"record group {group}: matrice non quadrata {matrix.shape}")
    return matrix


def relative_norm(value: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.linalg.norm(reference))
    return float(np.linalg.norm(value)) / max(denominator, np.finfo(float).tiny)


def load_matrix_cases(config: dict[str, Any], config_root: Path) -> list[MatrixCase]:
    model = config["model"]
    modes = list(model["selected_modes"])
    expected_shape = (len(modes), len(modes))
    expected_hash = model["baseline_fem_sha256"]
    expected_basis = model["basis_id"]
    tolerance = float(config["fit"]["symmetry_relative_tolerance"])
    cases: list[MatrixCase] = []
    units: str | None = None
    for entry in config["matrix_cases"]:
        directory = resolve(config_root, entry["directory"])
        matrix_path = directory / "kh_fixed_basis.csv"
        metadata_path = directory / "metadata.json"
        missing = [str(p) for p in (matrix_path, metadata_path) if not p.is_file()]
        if missing:
            raise FileNotFoundError(
                "mancano gli output Nastran proiettati:\n  " + "\n  ".join(missing)
                + "\nVedere matrices/README.md."
            )
        matrix = np.loadtxt(matrix_path, delimiter=",", comments="#", ndmin=2)
        metadata = json.loads(metadata_path.read_text())
        load_factor = float(entry["load_factor"])
        if matrix.shape != expected_shape or not np.isfinite(matrix).all():
            raise ValueError(f"{matrix_path}: attesa matrice finita {expected_shape}, trovata {matrix.shape}")
        checks = {
            "load_factor": math.isclose(float(metadata.get("load_factor", math.nan)), load_factor),
            "basis_id": metadata.get("basis_id") == expected_basis,
            "modes": metadata.get("modes") == modes,
            "fem_sha256": metadata.get("fem_sha256") == expected_hash,
            "projection": "Phi0" in str(metadata.get("projection", "")),
        }
        failed = [name for name, passed in checks.items() if not passed]
        if failed:
            raise ValueError(f"{metadata_path}: metadati incompatibili: {', '.join(failed)}")
        case_units = str(metadata.get("units", "")).strip()
        if not case_units:
            raise ValueError(f"{metadata_path}: units mancanti")
        if units is None:
            units = case_units
        elif units != case_units:
            raise ValueError(f"unita non uniformi: {units!r} contro {case_units!r}")
        symmetry_error = relative_norm(matrix - matrix.T, matrix)
        if symmetry_error > tolerance:
            raise ValueError(
                f"{matrix_path}: errore relativo di simmetria {symmetry_error:.3e} > {tolerance:.3e}"
            )
        cases.append(MatrixCase(load_factor, 0.5 * (matrix + matrix.T), metadata, directory, symmetry_error))
    cases.sort(key=lambda item: item.load_factor)
    if len({case.load_factor for case in cases}) != len(cases):
        raise ValueError("fattori di carico duplicati")
    return cases


def fit_slope(cases: list[MatrixCase], baseline_n: float) -> tuple[np.ndarray, MatrixCase, list[dict[str, float]]]:
    baseline = next((case for case in cases if math.isclose(case.load_factor, baseline_n)), None)
    if baseline is None:
        raise ValueError(f"manca il caso baseline n={baseline_n:g}")
    nonbaseline = [case for case in cases if not math.isclose(case.load_factor, baseline_n)]
    if not nonbaseline:
        raise ValueError("serve almeno un caso oltre il baseline")
    denominator = sum((case.load_factor - baseline_n) ** 2 for case in nonbaseline)
    slope = sum(
        ((case.load_factor - baseline_n) * (case.matrix - baseline.matrix) for case in nonbaseline),
        start=np.zeros_like(baseline.matrix),
    ) / denominator
    residuals = []
    for case in nonbaseline:
        prediction = baseline.matrix + (case.load_factor - baseline_n) * slope
        residuals.append({
            "load_factor": case.load_factor,
            "relative_residual": relative_norm(case.matrix - prediction, case.matrix),
        })
    return slope, baseline, residuals


def fmt(value: float) -> str:
    if abs(value) < 5e-16:
        value = 0.0
    return f"{value:.16e}"


def constants_text(config: dict[str, Any], slope: np.ndarray) -> str:
    modes = list(config["model"]["selected_modes"])
    mbdyn = config["mbdyn"]
    lines = [
        "# Generated by prestress_rom.py; do not edit.",
        "# KHN is d(K_h)/d(n_z) in the fixed Phi0 basis.",
        f'set: const integer {mbdyn["force_label"]} = {int(mbdyn["force_numeric_label"])};',
    ]
    for name, label in mbdyn["drive_numeric_labels"].items():
        lines.append(f"set: const integer {name} = {int(label)};")
    for i, mode_i in enumerate(modes):
        for j, mode_j in enumerate(modes):
            lines.append(f"set: const real KHN_{mode_i:02d}_{mode_j:02d} = {fmt(float(slope[i, j]))};")
    return "\n".join(lines) + "\n"


def force_text(config: dict[str, Any], slope: np.ndarray) -> str:
    model, mbdyn = config["model"], config["mbdyn"]
    modes = list(model["selected_modes"])
    qeq = list(mbdyn["q_equilibrium_constants"])
    drive = mbdyn["delta_n_drive"]
    rows = []
    for i, mode_i in enumerate(modes):
        terms = []
        for j, mode_j in enumerate(modes):
            # Keep zero terms too: the generated topology is invariant and auditable.
            terms.append(
                f'            element, {model["modal_joint"]}, joint, string, "q[{mode_j}]",',
            )
            terms[-1] += (
                f'\n                string, "-model::drive({drive},Time)*'
                f'KHN_{mode_i:02d}_{mode_j:02d}*(Var-{qeq[j]})"'
            )
        rows.append("        array, 6,\n" + ",\n".join(terms))
    return (
        "# Generated by prestress_rom.py; PRESTRESS_DELTA_N_DRIVE must already exist.\n"
        f'force: {mbdyn["force_label"]}, modal, {model["modal_joint"]},\n'
        f"    list, 6, {', '.join(str(mode) for mode in modes)},\n"
        + ",\n".join(rows)
        + ",\n    output, no;\n"
    )


def build(config_path: Path | None = None, output: Path | None = None) -> dict[str, Any]:
    config, config_root = load_config(config_path)
    model = config["model"]
    fem = resolve(config_root, model["baseline_fem"])
    actual_hash = sha256(fem)
    if actual_hash != model["baseline_fem_sha256"]:
        raise ValueError(f"hash FEM cambiato: {actual_hash}")
    cases = load_matrix_cases(config, config_root)
    slope, baseline, residuals = fit_slope(cases, float(model["baseline_load_factor"]))
    modes = list(model["selected_modes"])
    full_k = read_fem_matrix(fem, 10)
    indices = np.asarray([mode - 1 for mode in modes], dtype=int)
    fem_k = full_k[np.ix_(indices, indices)]
    baseline_mismatch = relative_norm(baseline.matrix - fem_k, baseline.matrix)
    effective_spectra = []
    for case in cases:
        effective = fem_k + (case.load_factor - float(model["baseline_load_factor"])) * slope
        eigvals = np.linalg.eigvalsh(0.5 * (effective + effective.T))
        effective_spectra.append({
            "load_factor": case.load_factor,
            "eigenvalues": eigvals.tolist(),
            "positive_definite": bool(np.all(eigvals > 0.0)),
        })
    output_dir = (output or config_root / "generated").expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "equation": "Q_prestress=-(n_z-1)*K_h_n*(q-q_eq)",
        "interpretation": "incremental_baseline_preserving",
        "basis_id": model["basis_id"],
        "modes": modes,
        "baseline_fem": str(fem),
        "baseline_fem_sha256": actual_hash,
        "matrix_units": cases[0].metadata["units"],
        "load_factor_definition": config["fit"]["load_factor_definition"],
        "load_factors": [case.load_factor for case in cases],
        "symmetry_errors": [case.symmetry_error for case in cases],
        "fit_residuals": residuals,
        "linearity_validated": len(cases) >= 3,
        "linearity_warning": (
            len(cases) < 3
            or max((item["relative_residual"] for item in residuals), default=0.0)
            > float(config["fit"]["linearity_warning_relative"])
        ),
        "baseline_kh_vs_fem_relative_mismatch": baseline_mismatch,
        "baseline_match_warning": baseline_mismatch > float(config["fit"]["baseline_match_warning_relative"]),
        "k_h_1g": baseline.matrix.tolist(),
        "k_h_n": slope.tolist(),
        "k_fem_active": fem_k.tolist(),
        "effective_stiffness_spectra": effective_spectra,
    }
    (output_dir / "prestress_rom.json").write_text(json.dumps(report, indent=2) + "\n")
    (output_dir / "prestress_constants.mbd").write_text(constants_text(config, slope))
    (output_dir / "prestress_modal_force.mbd").write_text(force_text(config, slope))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = build(args.config, args.output)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    print(json.dumps({
        "generated": str((args.output or ROOT / "generated").resolve()),
        "load_factors": report["load_factors"],
        "linearity_validated": report["linearity_validated"],
        "baseline_kh_vs_fem_relative_mismatch": report["baseline_kh_vs_fem_relative_mismatch"],
    }, indent=2))


if __name__ == "__main__":
    main()
