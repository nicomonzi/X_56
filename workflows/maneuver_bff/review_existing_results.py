#!/usr/bin/env python3
"""Read-only audit of source results; write separate diagnostic JSON only.

Does not prepare decks, rerun solvers, or change previous analysis artifacts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from analyse_prestressed_modes import parse_f06


ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_eigenvalues(root: Path) -> list[dict]:
    paths = sorted(root.glob("*/*.f06"))
    if not paths:
        raise RuntimeError(f"No F06 files under {root}")
    rows = []
    for path in paths:
        frequencies, _, eigenvalues = parse_f06(path, require_shapes=False)
        roots = [
            {
                "mode": mode,
                "eigenvalue_s_minus_2": value,
                "printed_cycles_hz": frequencies[mode],
                "negative": value < 0,
            }
            for mode, value in sorted(eigenvalues.items())
        ]
        rows.append({
            "source": str(path),
            "sha256": sha256(path),
            "roots": roots,
            "negative_roots_without_frequency_filter": [r for r in roots if r["negative"]],
            "interpretation": "Sign reported without suppression; physical cause not diagnosed.",
        })
    return rows


def audit_dlm(root: Path, margin: float) -> list[dict]:
    paths = sorted(root.glob("*.nc"))
    if not paths:
        raise RuntimeError(f"No NetCDF files under {root}")
    rows = []
    for path in paths:
        deck = path.with_suffix(".mbd")
        constants = {
            name: float(value) for name, value in re.findall(
                r"set:\s*const\s+real\s+(\w+)\s*=\s*([-+\d.eE]+)\s*;",
                deck.read_text(),
            )
        }
        k = constants["NASTRAN_DLM_ROM_K7"]
        c = constants["NASTRAN_DLM_ROM_C7"]
        qeq = constants["NASTRAN_DLM_ROM_Q7_EQ"]
        ramp_time = constants["NASTRAN_DLM_ROM_RAMP_TIME"]
        start = constants["SAS_OFF_START"] + margin
        end = constants["SAS_OFF_START"] + constants["SAS_OFF_DURATION"] - margin
        with Dataset(path) as data:
            time = np.asarray(data["time"][:])
            selected = (time >= start - 1e-8) & (time <= end + 1e-8)
            t = time[selected]
            if len(t) < 50:
                raise RuntimeError(f"Insufficient samples: {path}")
            q = np.asarray(data["elem.joint.5.a"][selected, 0])
            qd = np.asarray(data["elem.joint.5.aPrime"][selected, 0])
            actual = np.minimum(1.0, t / ramp_time) * (k * (q - qeq) + c * qd)
            legacy = 146.37892 * (q + 0.4304727364) + 3.67894 * qd
        delta = float(actual.mean() - legacy.mean())
        rows.append({
            "netcdf": str(path),
            "netcdf_sha256": sha256(path),
            "deck": str(deck),
            "deck_sha256": sha256(deck),
            "actual_coefficients": {"K7": k, "C7": c, "q7_eq": qeq, "ramp_time_s": ramp_time},
            "legacy_coefficients": {"K7": 146.37892, "C7": 3.67894, "q7_eq": -0.4304727364},
            "window_s": [float(t[0]), float(t[-1])],
            "samples": len(t),
            "legacy_Q7_mean": float(legacy.mean()),
            "actual_Q7_mean": float(actual.mean()),
            "actual_Q7_std": float(actual.std()),
            "actual_minus_legacy_Q7_mean": delta,
            "absolute_difference_over_actual_mean_percent": abs(delta / float(actual.mean())) * 100,
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nastran-root", type=Path, default=Path("/home/nicomonzi/ZENO/prestress_loads"))
    parser.add_argument("--recovery-cases", type=Path, default=Path("/mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/load_recovery/cases"))
    parser.add_argument("--margin-s", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=ROOT / "results/model_review/existing_results_audit.json")
    args = parser.parse_args()
    result = {
        "scope": "Diagnostic only; existing decks and previous artifacts unmodified.",
        "nastran": audit_eigenvalues(args.nastran_root),
        "dlm_load_recovery": audit_dlm(args.recovery_cases, args.margin_s),
        "decision": "Previous stable/robust-to-DLM gate not established; resolve audit findings first.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "negative_root_cases": [r["source"] for r in result["nastran"] if r["negative_roots_without_frequency_filter"]],
        "recovery_cases_checked": len(result["dlm_load_recovery"]),
    }, indent=2))


if __name__ == "__main__":
    main()
