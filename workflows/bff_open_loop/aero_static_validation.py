#!/usr/bin/env python3
"""Run rigid MBDyn X-56 polar points and compare them with NASTRAN.

The NASTRAN CMY printed in x56_polar.f06 is referred to the basic origin.
This script shifts it to the configured X-56 CG before any comparison.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/bff_aero_validation_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

import analyse_open_loop as analysis


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[1]
MODEL = ROOT / "main_bff_open_loop.mbd"
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")
NASTRAN_JSON = REPO_ROOT / "validation/aero_polar/test/nastran_coefficients.json"
X_CG_IN = 163.187383385809
Y_CG_IN = 0.110529571088
Z_CG_IN = 101.239797358848
S_REF_IN2 = 8064.0
C_REF_IN = 24.0
RHO_SLUG_IN3 = 9.7284e-8


def replace_constant(text: str, name: str, expression: str) -> str:
    pattern = rf"(?m)^(\s*set:\s*const\s+real\s+{re.escape(name)}\s*=\s*)[^;]+;"
    value, count = re.subn(pattern, rf"\g<1>{expression};", text)
    if count != 1:
        raise RuntimeError(f"constant {name}: expected one definition, found {count}")
    return value


def prepare_case(angle_deg: float, velocity_mps: float, body_surface_deg: float, wing_surface_deg: float, case_dir: Path) -> tuple[Path, Path]:
    case_dir.mkdir(parents=True, exist_ok=True)
    include_dir = case_dir / "INCLUDE"
    if include_dir.exists():
        shutil.rmtree(include_dir)
    shutil.copytree(ROOT / "INCLUDE", include_dir)
    # A rigid static derivative must test the steady C81 kernel.  The
    # Theodorsen wrapper has the same t->infinity value but needs time to
    # converge, while the retained modal coordinates would then acquire a
    # static aeroelastic deformation and contaminate the rigid comparison.
    aerobody = include_dir / "aerobody.mbd"
    aerobody.write_text(aerobody.read_text().replace("theodorsen, c81", "c81"))

    constants = (include_dir / "setconst.mbd").read_text()
    constants = replace_constant(constants, "TRIM_PITCH", f"({angle_deg:.9f})*deg2rad")
    constants = replace_constant(constants, "TRIM_SURFACE", f"({wing_surface_deg:.9f})*deg2rad")
    constants = replace_constant(constants, "LONG_TRIM_BIAS", "0.*deg2rad")
    constants = replace_constant(constants, "BODY_TRIM_CORRECTION", "0.*deg2rad")
    constants = replace_constant(constants, "BODY_TRIM_SURFACE", f"({body_surface_deg:.9f})*deg2rad")
    for mode in range(7, 13):
        constants = replace_constant(constants, f"Q_INIT_{mode}", "0.")
    (include_dir / "setconst.mbd").write_text(constants)

    text = MODEL.read_text()
    text = replace_constant(text, "V_INF", f"{velocity_mps:.9f}")
    text = replace_constant(text, "FINAL_TIME", "0.05")
    text = replace_constant(text, "TIME_STEP", "0.01")
    text = text.replace('"./INCLUDE/', f'"{include_dir}/')
    text = text.replace(
        "rotation orientation, reference, global, eye,\n"
        "        position, reference, global, X_CG, Y_CG, Z_CG,\n"
        "        position orientation, reference, global, eye,\n"
        "        rotation orientation, reference, global, eye,\n"
        "        position constraint, active, active, inactive, null,\n"
        "        orientation constraint, inactive, inactive, inactive, null;",
        "rotation orientation, reference, global, euler123, 0., TRIM_PITCH, 0.,\n"
        "        position, reference, global, X_CG, Y_CG, Z_CG,\n"
        "        position orientation, reference, global, eye,\n"
        "        rotation orientation, reference, global, euler123, 0., TRIM_PITCH, 0.,\n"
        "        position constraint, active, active, active, null,\n"
        "        orientation constraint, active, active, active, null;",
    )
    if "orientation constraint, active, active, active" not in text:
        raise RuntimeError("failed to make the rigid polar constraint")
    stem = f"alpha_{angle_deg:+06.2f}_body_{body_surface_deg:+06.2f}_wing_{wing_surface_deg:+06.2f}".replace("+", "p").replace("-", "m").replace(".", "p")
    input_path = case_dir / f"{stem}.mbd"
    input_path.write_text(text)
    return input_path, case_dir / stem


def integrate_aerodynamics(nc_path: Path, velocity_mps: float) -> dict[str, float]:
    cg = np.array([X_CG_IN, Y_CG_IN, Z_CG_IN])
    force = np.zeros(3)
    moment = np.zeros(3)
    with Dataset(nc_path) as data:
        # Use the initial output before the flexible modal coordinates react.
        # The copied polar deck uses the steady C81 kernel (the t->infinity
        # limit of the Theodorsen wrapper), so this is the rigid derivative.
        index = 0
        for label, span in analysis.aerodynamic_spans().items():
            for point, weight in enumerate(analysis.GAUSS_W):
                scale = 0.5 * span * weight
                f = np.asarray(data[f"elem.aerodynamic.{label}.F_{point}"][index]) * scale
                m = np.asarray(data[f"elem.aerodynamic.{label}.M_{point}"][index]) * scale
                x = np.asarray(data[f"elem.aerodynamic.{label}.X_{point}"][index])
                force += f
                moment += m + np.cross(x - cg, f)
    speed_inps = velocity_mps / analysis.IN_TO_M
    q_dyn = 0.5 * RHO_SLUG_IN3 * speed_inps**2
    return {
        "CX": float(force[0] / (q_dyn * S_REF_IN2)),
        "CY": float(force[1] / (q_dyn * S_REF_IN2)),
        "CZ": float(force[2] / (q_dyn * S_REF_IN2)),
        "CMX_CG": float(moment[0] / (q_dyn * S_REF_IN2 * C_REF_IN)),
        "CMY_CG": float(moment[1] / (q_dyn * S_REF_IN2 * C_REF_IN)),
        "CMZ_CG": float(moment[2] / (q_dyn * S_REF_IN2 * C_REF_IN)),
    }


def nastran_polar() -> dict[float, dict[str, float]]:
    source = json.loads(NASTRAN_JSON.read_text())
    result = {}
    for key, row in source.items():
        angle = round(float(key))
        # M_CG = M_origin - r_CG x F.  In the x-z plane this gives
        # CMY_CG = CMY_origin + x_CG/c*CZ - z_CG/c*CX.
        corrected = row["CMY"] + (X_CG_IN / C_REF_IN) * row["CZ"] - (Z_CG_IN / C_REF_IN) * row["CX"]
        result[angle] = {
            "CZ": float(row["CZ"]),
            "CX": float(row["CX"]),
            "CMY_CG": float(corrected),
            "CMY_origin": float(row["CMY"]),
        }
    return result


def fit_slope(rows: dict[float, dict[str, float]], key: str) -> float | None:
    angles = np.array(sorted(a for a in rows if abs(a) <= 2.0))
    if len(angles) < 2 or np.ptp(angles) == 0.0:
        return None
    values = np.array([rows[float(a)][key] for a in angles])
    return float(np.polyfit(angles, values, 1)[0])


def plot_results(mbdyn: dict[float, dict[str, float]], nastran: dict[float, dict[str, float]], output: Path) -> None:
    angles = np.array(sorted(mbdyn))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for axis, key, ylabel in zip(axes, ("CZ", "CX", "CMY_CG"), ("CZ", "CX", "CMY about CG")):
        axis.plot(angles, [mbdyn[a][key] for a in angles], "o-", label="MBDyn strip/C81")
        axis.plot(angles, [nastran[a][key] for a in angles], "s--", label="NASTRAN DLM")
        axis.axhline(0.0, color="0.3", lw=0.7)
        axis.set_xlabel("angle of attack [deg]")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
    axes[0].legend()
    fig.savefig(output / "nastran_mbdyn_rigid_polar.png", dpi=190)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--velocity", type=float, default=60.8421)
    parser.add_argument("--angles", type=float, nargs="+", default=[-2.0, -1.0, 0.0, 1.0, 2.0])
    parser.add_argument("--body-surface", type=float, default=0.0)
    parser.add_argument("--wing-surface", type=float, default=0.0)
    parser.add_argument("--output", type=Path, default=ROOT / "aero_validation")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    args.output = args.output.resolve()
    case_dir = args.output / "mbdyn_cases"
    mbdyn_rows: dict[float, dict[str, float]] = {}
    for angle in args.angles:
        input_path, prefix = prepare_case(angle, args.velocity, args.body_surface, args.wing_surface, case_dir)
        nc_path = prefix.with_suffix(".nc")
        if args.overwrite or not nc_path.exists():
            result = subprocess.run(
                [str(MBDYN), "-s", "-f", str(input_path), "-o", str(prefix)],
                cwd=ROOT, text=True, capture_output=True,
            )
            prefix.with_suffix(".stdout").write_text(result.stdout + result.stderr)
            if result.returncode or not nc_path.exists():
                raise RuntimeError(f"MBDyn failed for alpha={angle}; see {prefix.with_suffix('.stdout')}")
        mbdyn_rows[float(angle)] = integrate_aerodynamics(nc_path, args.velocity)
    nastran = nastran_polar()
    selected_nastran = {float(a): nastran[round(a)] for a in args.angles}
    payload = {
        "reference": {"S_ref_in2": S_REF_IN2, "c_ref_in": C_REF_IN, "CG_in": [X_CG_IN, Y_CG_IN, Z_CG_IN]},
        "surface_settings_deg": {"body": args.body_surface, "wing": args.wing_surface},
        "warning": "NASTRAN CMY shifted from basic origin to CG; legacy synthetic comparison is invalid",
        "mbdyn": mbdyn_rows,
        "nastran": selected_nastran,
        "slopes_per_deg": {
            "mbdyn_CZ_alpha": fit_slope(mbdyn_rows, "CZ"),
            "nastran_CZ_alpha": fit_slope(selected_nastran, "CZ"),
            "mbdyn_CMY_CG_alpha": fit_slope(mbdyn_rows, "CMY_CG"),
            "nastran_CMY_CG_alpha": fit_slope(selected_nastran, "CMY_CG"),
        },
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "aero_comparison.json").write_text(json.dumps(payload, indent=2))
    plot_results(mbdyn_rows, selected_nastran, args.output)
    print(json.dumps(payload["slopes_per_deg"], indent=2))


if __name__ == "__main__":
    main()
