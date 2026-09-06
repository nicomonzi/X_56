#!/usr/bin/env python3
"""Average MBDyn aerodynamic loads in the aircraft frame and write Nastran cards."""

from __future__ import annotations

import argparse
import csv
import functools
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from campaign import BASELINE_DIR, CONFIG


DEFAULT_RESULTS = Path(CONFIG["execution"]["output_root"]) / "load_recovery"
AEROBODY = BASELINE_DIR / "INCLUDE/aerobody.mbd"
CONTROL_PARENT = {
    880004: 990004,
    880008: 990008,
    880011: 990011,
    880014: 990014,
    880017: 990017,
    881004: 991004,
    881008: 991008,
    881011: 991011,
    881014: 991014,
    881017: 991017,
}
GAUSS_XI = math.sqrt(3.0 / 5.0)
GAUSS_WEIGHTS = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
G_IN_S2 = float(CONFIG["physics"]["gravity_m_s2"]) * 100.0 / 2.54
WTMASS = 0.002591
MASS_LBM = 419.1477
CG_NASTRAN_IN = np.array([163.1974, 0.1184380, 101.2380])
SUPPORT_GRID = 10062
SUPPORT_POSITION_IN = np.array([163.170, 0.0, 100.445])
INERTIA_CG_LBM_IN2 = np.array([
    [1.616757e6, -2.007283e3, 1.302301e4],
    [-2.007283e3, 2.383475e5, -7.370287e2],
    [1.302301e4, -7.370287e2, 1.831537e6],
])
NASTRAN_SOURCE = Path("/home/nicomonzi/X_56/validation/modal/sol103_60_modes")
BASELINE_FEM = Path("/home/nicomonzi/X_56/workflows/bff_open_loop/INCLUDE/mbdyn_modal.fem")


def _constant(text: str, name: str) -> float:
    match = re.search(
        rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*([-+0-9.eE]+)\s*;",
        text,
    )
    if not match:
        raise RuntimeError(f"costante {name} non trovata")
    return float(match.group(1))


def _metadata(text: str) -> dict:
    match = re.search(r"(?m)^# STIFFNESS_STUDY_METADATA (\{.*\})$", text)
    if not match:
        raise RuntimeError("metadata del caso non trovati")
    return json.loads(match.group(1))


def element_nodes() -> dict[int, int]:
    pairs = re.findall(
        r"(?m)^\s*aerodynamic body:\s*(\d+)\s*,\s*(\d+)\s*,",
        AEROBODY.read_text(),
    )
    result = {int(element): int(node) for element, node in pairs}
    if set(result) != set(range(1, 59)):
        raise RuntimeError("mappa aerodinamica incompleta o duplicata")
    return result


@functools.lru_cache(maxsize=1)
def canonical_mode7_load() -> dict:
    """Return the unit generalized load M_lumped phi_7 on every FEM grid."""
    labels: list[int] = []
    coordinates = {5: [], 6: [], 7: []}
    mode_shape: list[list[float]] = []
    mass_diagonal: list[list[float]] = []
    group = 0
    mode: int | None = None
    with BASELINE_FEM.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            marker = re.search(r"RECORD GROUP\s+(\d+)", line)
            if marker:
                group = int(marker.group(1))
                continue
            mode_marker = re.search(r"NORMAL MODE SHAPE\s*#\s*(\d+)", line)
            if mode_marker:
                mode = int(mode_marker.group(1))
                continue
            if line.lstrip().startswith("*") or not line.split():
                continue
            fields = line.split()
            if group == 2:
                labels.extend(int(item) for item in fields)
            elif group in coordinates:
                coordinates[group].extend(float(item) for item in fields)
            elif group == 8 and mode == 7 and len(fields) == 6:
                mode_shape.append([float(item) for item in fields])
            elif group == 11 and len(fields) == 6:
                mass_diagonal.append([float(item) for item in fields])
    phi = np.asarray(mode_shape)
    mdiag = np.asarray(mass_diagonal)
    xyz = np.column_stack([coordinates[index] for index in (5, 6, 7)])
    expected = (len(labels), 6)
    if phi.shape != expected or mdiag.shape != expected or xyz.shape != (len(labels), 3):
        raise RuntimeError(
            f"FEM incompleto per la distribuzione DLM: labels={len(labels)}, "
            f"phi={phi.shape}, mdiag={mdiag.shape}, xyz={xyz.shape}"
        )
    modal_projection = float(np.sum(phi * mdiag * phi))
    if not 0.98 <= modal_projection <= 1.02:
        raise RuntimeError(f"normalizzazione modo 7 inattesa: {modal_projection}")
    unit_load = mdiag * phi / modal_projection
    force = unit_load[:, :3]
    moment = unit_load[:, 3:]
    resultant_force = np.sum(force, axis=0)
    resultant_moment_cg = np.sum(
        moment + np.cross(xyz - CG_NASTRAN_IN, force), axis=0
    )
    return {
        "labels": np.asarray(labels, dtype=int),
        "coordinates": xyz,
        "unit_load": unit_load,
        "raw_modal_projection": modal_projection,
        "normalized_modal_projection": float(np.sum(phi * unit_load)),
        "unit_resultant_force": resultant_force,
        "unit_resultant_moment_about_cg": resultant_moment_cg,
    }


def rotation_matrices(rotvec: np.ndarray) -> np.ndarray:
    """Rodrigues matrices for MBDyn orientation vectors."""
    result = np.empty((len(rotvec), 3, 3))
    identity = np.eye(3)
    for index, vector in enumerate(rotvec):
        angle = float(np.linalg.norm(vector))
        if angle < 1e-14:
            result[index] = identity
            continue
        axis = vector / angle
        skew = np.array([
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ])
        result[index] = identity + math.sin(angle) * skew + (1.0 - math.cos(angle)) * (skew @ skew)
    return result


def to_body(rotation: np.ndarray, vector: np.ndarray) -> np.ndarray:
    return np.einsum("tji,tj->ti", rotation, vector)


def extract_case(nc_path: Path, margin_s: float) -> dict:
    mbd_path = nc_path.with_suffix(".mbd")
    text = mbd_path.read_text()
    metadata = _metadata(text)
    start = _constant(text, "SAS_OFF_START") + margin_s
    end = _constant(text, "SAS_OFF_START") + _constant(text, "SAS_OFF_DURATION") - margin_s
    mapping = element_nodes()

    with Dataset(nc_path) as data:
        time = np.asarray(data["time"][:], dtype=float)
        mask = (time >= start - 1e-8) & (time <= end + 1e-8)
        if np.count_nonzero(mask) < 50:
            raise RuntimeError(f"campioni insufficienti nella finestra {start:.3f}-{end:.3f} s")
        t = time[mask]
        rotation = rotation_matrices(np.asarray(data["node.struct.990000.Phi"][mask], dtype=float))
        base_x = np.asarray(data["node.struct.990000.X"][mask], dtype=float)
        base_omega = np.asarray(data["node.struct.990000.Omega"][mask], dtype=float)
        base_velocity = np.asarray(data["node.struct.990000.XP"][mask], dtype=float)
        # Modal interface nodes do not expose XPP/OmegaP in this MBDyn build.
        # The fixed 50 Hz output sampling is sufficient for the quasi-steady
        # mean needed here; second-order finite differences avoid another
        # model approximation and the standard deviations quantify ripple.
        base_xpp = np.gradient(base_velocity, t, axis=0, edge_order=2)
        base_omegap = np.gradient(base_omega, t, axis=0, edge_order=2)

        gravity_global = np.zeros_like(base_xpp)
        gravity_global[:, 2] = -G_IN_S2
        specific_g_body = to_body(rotation, base_xpp - gravity_global) / G_IN_S2
        angular_accel_body = to_body(rotation, base_omegap)
        angular_velocity_body = to_body(rotation, base_omega)

        target_positions = {}
        for source in set(mapping.values()):
            target = CONTROL_PARENT.get(source, source)
            if target not in target_positions:
                target_positions[target] = np.asarray(
                    data[f"node.struct.{target}.X"][mask], dtype=float
                )

        node_forces = {node: np.zeros((len(t), 3)) for node in target_positions}
        node_moments = {node: np.zeros((len(t), 3)) for node in target_positions}
        total_force = np.zeros((len(t), 3))
        total_moment_cg = np.zeros((len(t), 3))

        for element, source_node in mapping.items():
            target = CONTROL_PARENT.get(source_node, source_node)
            x_gp = [np.asarray(data[f"elem.aerodynamic.{element}.X_{gp}"][mask], dtype=float) for gp in range(3)]
            f_gp = [np.asarray(data[f"elem.aerodynamic.{element}.F_{gp}"][mask], dtype=float) for gp in range(3)]
            m_gp = [np.asarray(data[f"elem.aerodynamic.{element}.M_{gp}"][mask], dtype=float) for gp in range(3)]
            half_span = np.linalg.norm(x_gp[2] - x_gp[0], axis=1) / (2.0 * GAUSS_XI)
            for gp, weight in enumerate(GAUSS_WEIGHTS):
                scale = (weight * half_span)[:, None]
                force_body = to_body(rotation, f_gp[gp])
                moment_body = to_body(rotation, m_gp[gp])
                arm_node = to_body(rotation, x_gp[gp] - target_positions[target])
                arm_cg = to_body(rotation, x_gp[gp] - base_x)
                integrated_force = scale * force_body
                node_forces[target] += integrated_force
                node_moments[target] += scale * (moment_body + np.cross(arm_node, force_body))
                total_force += integrated_force
                total_moment_cg += scale * (moment_body + np.cross(arm_cg, force_body))

        rows = []
        for node in sorted(node_forces):
            f_mean = np.mean(node_forces[node], axis=0)
            m_mean = np.mean(node_moments[node], axis=0)
            f_std = np.std(node_forces[node], axis=0)
            m_std = np.std(node_moments[node], axis=0)
            rows.append({
                "grid": node,
                "fx_lbf": f_mean[0], "fy_lbf": f_mean[1], "fz_lbf": f_mean[2],
                "mx_lbf_in": m_mean[0], "my_lbf_in": m_mean[1], "mz_lbf_in": m_mean[2],
                "fx_std_lbf": f_std[0], "fy_std_lbf": f_std[1], "fz_std_lbf": f_std[2],
                "mx_std_lbf_in": m_std[0], "my_std_lbf_in": m_std[1], "mz_std_lbf_in": m_std[2],
            })

        q = np.asarray(data["elem.joint.5.a"][mask, 0], dtype=float)
        qp = np.asarray(data["elem.joint.5.aPrime"][mask, 0], dtype=float)
        k_dlm = 146.37892000
        c_dlm = 3.67894000
        q_eq = -0.4304727364
        dlm_q7 = k_dlm * (q - q_eq) + c_dlm * qp

    return {
        "source_netcdf": str(nc_path),
        "metadata": metadata,
        "window_s": [float(t[0]), float(t[-1])],
        "sample_count": len(t),
        "specific_acceleration_body_g_mean": np.mean(specific_g_body, axis=0).tolist(),
        "specific_acceleration_body_g_std": np.std(specific_g_body, axis=0).tolist(),
        "angular_acceleration_body_rad_s2_mean": np.mean(angular_accel_body, axis=0).tolist(),
        "angular_acceleration_body_rad_s2_std": np.std(angular_accel_body, axis=0).tolist(),
        "angular_velocity_body_rad_s_mean": np.mean(angular_velocity_body, axis=0).tolist(),
        "angular_velocity_body_rad_s_std": np.std(angular_velocity_body, axis=0).tolist(),
        "aerodynamic_force_body_lbf_mean": np.mean(total_force, axis=0).tolist(),
        "aerodynamic_force_body_lbf_std": np.std(total_force, axis=0).tolist(),
        "aerodynamic_moment_about_cg_body_lbf_in_mean": np.mean(total_moment_cg, axis=0).tolist(),
        "aerodynamic_moment_about_cg_body_lbf_in_std": np.std(total_moment_cg, axis=0).tolist(),
        "dlm_mode7_generalized_force_mean": float(np.mean(dlm_q7)),
        "dlm_mode7_generalized_force_std": float(np.std(dlm_q7)),
        "dlm_distribution_note": (
            "The DLM correction is available only as generalized mode-7 force; "
            "it has no unique physical nodal load distribution. Separate Nastran decks exclude it "
            "or include the documented canonical M_lumped*phi_7 distribution."
        ),
        "nodal_loads": rows,
    }


def write_case(result: dict, output: Path) -> None:
    nominal_n = float(result["metadata"]["nominal_load_factor"])
    tag = f"n{nominal_n:.1f}".replace(".", "p")
    case_dir = output / tag
    case_dir.mkdir(parents=True, exist_ok=True)
    fields = list(result["nodal_loads"][0])
    with (case_dir / "nodal_loads.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(result["nodal_loads"])

    aero_sid = 610
    inertia_sid = 611
    balance_sid = 612
    total_sid = 613
    dlm_sid = 614
    total_with_dlm_sid = 615
    dlm_balance_sid = 616
    rotational_inertia_sid = 617
    combined_inertia_sid = 618
    specific_g = np.asarray(result["specific_acceleration_body_g_mean"], dtype=float)
    effective_accel = -specific_g * G_IN_S2
    omega = np.asarray(result["angular_velocity_body_rad_s_mean"], dtype=float)
    alpha = np.asarray(result["angular_acceleration_body_rad_s2_mean"], dtype=float)
    aero_force = np.asarray(result["aerodynamic_force_body_lbf_mean"], dtype=float)
    aero_moment_cg = np.asarray(
        result["aerodynamic_moment_about_cg_body_lbf_in_mean"], dtype=float
    )
    inertia_force = MASS_LBM * WTMASS * effective_accel
    inertia_moment = -WTMASS * (
        INERTIA_CG_LBM_IN2 @ alpha
        + np.cross(omega, INERTIA_CG_LBM_IN2 @ omega)
    )
    raw_force_residual = aero_force + inertia_force
    raw_moment_residual_cg = aero_moment_cg + inertia_moment
    balance_force = -raw_force_residual
    arm_support_cg = SUPPORT_POSITION_IN - CG_NASTRAN_IN
    balance_moment_support = (
        -raw_moment_residual_cg - np.cross(arm_support_cg, balance_force)
    )
    force_scale = sum(
        float(np.linalg.norm([row["fx_lbf"], row["fy_lbf"], row["fz_lbf"]]))
        for row in result["nodal_loads"]
    ) + float(np.linalg.norm(inertia_force))
    moment_scale = sum(
        float(np.linalg.norm([row["mx_lbf_in"], row["my_lbf_in"], row["mz_lbf_in"]]))
        for row in result["nodal_loads"]
    ) + float(np.linalg.norm(aero_moment_cg)) + float(np.linalg.norm(inertia_moment))
    result["nastran_balance"] = {
        "wtmass": WTMASS,
        "mass_lbm": MASS_LBM,
        "cg_basic_in": CG_NASTRAN_IN.tolist(),
        "support_grid": SUPPORT_GRID,
        "effective_gravity_acceleration_body_in_s2": effective_accel.tolist(),
        "distributed_inertia_force_lbf": inertia_force.tolist(),
        "rotational_inertia_moment_about_cg_lbf_in": inertia_moment.tolist(),
        "rforce_axis_grid_offset_from_cg_in": (SUPPORT_POSITION_IN - CG_NASTRAN_IN).tolist(),
        "rforce_angular_velocity_rad_s": omega.tolist(),
        "rforce_opposite_angular_acceleration_rad_s2": (-alpha).tolist(),
        "raw_force_residual_lbf": raw_force_residual.tolist(),
        "raw_moment_residual_about_cg_lbf_in": raw_moment_residual_cg.tolist(),
        "raw_force_residual_ratio": float(np.linalg.norm(raw_force_residual) / max(force_scale, 1e-30)),
        "raw_moment_residual_ratio": float(np.linalg.norm(raw_moment_residual_cg) / max(moment_scale, 1e-30)),
        "balancing_force_at_grid_10062_lbf": balance_force.tolist(),
        "balancing_moment_at_grid_10062_lbf_in": balance_moment_support.tolist(),
        "interpretation": (
            "The balancing wrench represents the small omitted rigid-body trim reactions "
            "(including the MBDyn x/y pin/thrust balance). It makes the static load set self-equilibrated."
        ),
    }
    canonical = canonical_mode7_load()
    dlm_q7 = float(result["dlm_mode7_generalized_force_mean"])
    dlm_load = canonical["unit_load"] * dlm_q7
    dlm_force_residual = canonical["unit_resultant_force"] * dlm_q7
    dlm_moment_residual_cg = canonical["unit_resultant_moment_about_cg"] * dlm_q7
    dlm_balance_force = -dlm_force_residual
    dlm_balance_moment_support = (
        -dlm_moment_residual_cg
        - np.cross(arm_support_cg, dlm_balance_force)
    )
    result["dlm_canonical_distribution"] = {
        "definition": "f = M_lumped * phi_7 * Q_7 / (phi_7^T M_lumped phi_7)",
        "raw_modal_projection": canonical["raw_modal_projection"],
        "normalized_modal_projection": canonical["normalized_modal_projection"],
        "generalized_force_Q7": dlm_q7,
        "force_resultant_before_balance_lbf": dlm_force_residual.tolist(),
        "moment_resultant_about_cg_before_balance_lbf_in": dlm_moment_residual_cg.tolist(),
        "limitation": (
            "This is a canonical modal-equivalent distribution, not a unique DLM panel load. "
            "Use the no-DLM and with-DLM decks as an uncertainty bracket."
        ),
    }
    (case_dir / "load_summary.json").write_text(json.dumps(result, indent=2))

    cards = [
        "$ Self-equilibrated MBDyn preload, averaged in the aircraft/body frame.",
        f"$ Source: {result['source_netcdf']}",
        f"$ Window: {result['window_s'][0]:.6f} to {result['window_s'][1]:.6f} s",
        "$ Units: lbf, inch, second. Loads act at the 45 RBE3 reference grids.",
        "$ GRAV is the measured gravity-minus-translation-acceleration field.",
    ]
    for row in result["nodal_loads"]:
        grid = int(row["grid"])
        force = [row["fx_lbf"], row["fy_lbf"], row["fz_lbf"]]
        moment = [row["mx_lbf_in"], row["my_lbf_in"], row["mz_lbf_in"]]
        cards.append(
            f"FORCE,{aero_sid},{grid},0,1.0,{force[0]:.10e},{force[1]:.10e},{force[2]:.10e}"
        )
        cards.append(
            f"MOMENT,{aero_sid},{grid},0,1.0,{moment[0]:.10e},{moment[1]:.10e},{moment[2]:.10e}"
        )
    cards.extend([
        f"GRAV,{inertia_sid},0,1.0,{effective_accel[0]:.10e},{effective_accel[1]:.10e},{effective_accel[2]:.10e}",
    ])
    omega_norm = float(np.linalg.norm(omega))
    alpha_norm = float(np.linalg.norm(alpha))
    if omega_norm > 1e-14:
        axis = omega / omega_norm
        cards.append(
            f"RFORCE,{rotational_inertia_sid},{SUPPORT_GRID},0,{omega_norm/(2.0*math.pi):.6e},"
            f"{axis[0]:.6e},{axis[1]:.6e},{axis[2]:.6e},1"
        )
    if alpha_norm > 1e-14:
        # RFORCE's RACC term acts with the angular acceleration; use the
        # opposite axis for the d'Alembert inertia load.
        axis = -alpha / alpha_norm
        cards.append(
            f"RFORCE,{rotational_inertia_sid},{SUPPORT_GRID},0,0.0,"
            f"{axis[0]:.6e},{axis[1]:.6e},{axis[2]:.6e},1,+RF617A"
        )
        cards.append(f"+RF617A,{alpha_norm/(2.0*math.pi):.6e}")
    cards.extend([
        f"LOAD,{combined_inertia_sid},1.0,1.0,{inertia_sid},1.0,{rotational_inertia_sid}",
        f"FORCE,{balance_sid},{SUPPORT_GRID},0,1.0,{balance_force[0]:.10e},{balance_force[1]:.10e},{balance_force[2]:.10e}",
        f"MOMENT,{balance_sid},{SUPPORT_GRID},0,1.0,{balance_moment_support[0]:.10e},{balance_moment_support[1]:.10e},{balance_moment_support[2]:.10e}",
        f"LOAD,{total_sid},1.0,1.0,{aero_sid},1.0,{combined_inertia_sid},1.0,{balance_sid}",
    ])
    cards.append("$ Canonical physical realization of the MBDyn generalized DLM mode-7 load.")
    for label, load in zip(canonical["labels"], dlm_load):
        # FORCE and MOMENT are independent Nastran cards.  A modal vector can
        # have only translations or only rotations at a grid; never emit a
        # card whose direction vector is exactly zero (UFM 9994).
        if float(np.linalg.norm(load[:3])) > 1e-15:
            cards.append(
                f"FORCE,{dlm_sid},{int(label)},0,1.0,{load[0]:.10e},{load[1]:.10e},{load[2]:.10e}"
            )
        if float(np.linalg.norm(load[3:])) > 1e-15:
            cards.append(
                f"MOMENT,{dlm_sid},{int(label)},0,1.0,{load[3]:.10e},{load[4]:.10e},{load[5]:.10e}"
            )
    cards.extend([
        f"FORCE,{dlm_balance_sid},{SUPPORT_GRID},0,1.0,{dlm_balance_force[0]:.10e},{dlm_balance_force[1]:.10e},{dlm_balance_force[2]:.10e}",
        f"MOMENT,{dlm_balance_sid},{SUPPORT_GRID},0,1.0,{dlm_balance_moment_support[0]:.10e},{dlm_balance_moment_support[1]:.10e},{dlm_balance_moment_support[2]:.10e}",
        f"LOAD,{total_with_dlm_sid},1.0,1.0,{total_sid},1.0,{dlm_sid},1.0,{dlm_balance_sid}",
    ])
    (case_dir / "preload_loads.bdf").write_text("\n".join(cards) + "\n")
    write_prestressed_modes_deck(case_dir, nominal_n, total_sid, "nodlm")
    write_prestressed_modes_deck(case_dir, nominal_n, total_with_dlm_sid, "with_dlm")


def write_prestressed_modes_deck(
    case_dir: Path, nominal_n: float, load_sid: int, variant: str
) -> None:
    """Create a self-contained two-subcase SOL103 following the MSC example."""
    bulk_source = NASTRAN_SOURCE / "BULK"
    bulk_target = case_dir / "BULK"
    if not bulk_target.exists():
        shutil.copytree(bulk_source, bulk_target)
    shutil.copy2(NASTRAN_SOURCE / "MAIN/rbe3s.bdf", case_dir / "rbe3s.bdf")
    source = (NASTRAN_SOURCE / "MAIN/sol103_60_modes.bdf").read_text()
    bulk = source.split("BEGIN BULK", 1)[1]
    bulk = re.sub(r"(?m)^\s*EIGRL\s+1.*$", "", bulk, count=1)
    bulk = bulk.replace("../BULK/", "BULK/")
    bulk = bulk.rsplit("ENDDATA", 1)[0].rstrip()
    header = f"""NASTRAN Q4TAPER=0.7,T3SKEW=3.0
INIT MASTER(S)
SOL 103
CEND
TITLE=PRESTIFF MODES V66.75 NNOM{nominal_n:.1f} {variant.upper()}
ECHO=NONE
SET 901 = 990001 THRU 990023, 991002 THRU 991023
SUBCASE 1
    SUBTITLE=MBDYN BODY-FRAME STATIC PRELOAD
    SPC=620
    LOAD={load_sid}
    DISPLACEMENT(PLOT,SORT1,REAL)=ALL
    SPCFORCES(PLOT,SORT1,REAL)=ALL
    STRESS(PLOT,SORT1,REAL,VONMISES,BILIN)=ALL
SUBCASE 2
    SUBTITLE=CG-SUPPORTED NORMAL MODES ABOUT PRELOADED STATE
    STATSUB=1
    SPC=620
    METHOD=1
    VECTOR(SORT1,REAL,PRINT,PLOT)=901
BEGIN BULK

EIGRL,1,,,30,0,,,MASS
PARAM,TESTNEG,3
SPC1,620,123456,{SUPPORT_GRID}
INCLUDE 'preload_loads.bdf'
"""
    (case_dir / f"prestressed_modes_{variant}.bdf").write_text(header + bulk + "\nENDDATA\n")
    (case_dir / "run_nastran.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        "solver=\"${NASTRAN_CMD:-nastran}\"\n"
        "\"$solver\" prestressed_modes_nodlm.bdf scr=yes old=no\n"
        "\"$solver\" prestressed_modes_with_dlm.bdf scr=yes old=no\n"
    )
    (case_dir / "run_nastran.sh").chmod(0o755)
    (case_dir / "run_zeno.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        "nast old=no prestressed_modes_nodlm.bdf\n"
        "nast old=no prestressed_modes_with_dlm.bdf\n"
    )
    (case_dir / "run_zeno.sh").chmod(0o755)
    if abs(nominal_n - 1.0) < 1e-12 and variant == "nodlm":
        unloaded_header = f"""NASTRAN Q4TAPER=0.7,T3SKEW=3.0
INIT MASTER(S)
SOL 103
CEND
TITLE=CG-SUPPORTED UNLOADED REFERENCE MODES
ECHO=NONE
SET 901 = 990001 THRU 990023, 991002 THRU 991023
SUBCASE 1
    SUBTITLE=CG-SUPPORTED UNLOADED NORMAL MODES
    SPC=620
    METHOD=1
    VECTOR(SORT1,REAL,PRINT,PLOT)=901
BEGIN BULK

EIGRL,1,,,30,0,,,MASS
SPC1,620,123456,{SUPPORT_GRID}
"""
        (case_dir / "supported_unloaded.bdf").write_text(
            unloaded_header + bulk + "\nENDDATA\n"
        )
        (case_dir / "run_supported_unloaded_zeno.sh").write_text(
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            "nast old=no supported_unloaded.bdf\n"
        )
        (case_dir / "run_supported_unloaded_zeno.sh").chmod(0o755)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, nargs="?", default=DEFAULT_RESULTS)
    parser.add_argument("--margin", type=float, default=0.25)
    args = parser.parse_args()
    if not 0.0 <= args.margin < 0.75:
        raise SystemExit("--margin deve essere compreso fra 0 e 0.75 s")
    results_dir = args.results.expanduser().resolve()
    output = results_dir / "prestress_loads"
    nc_files = sorted((results_dir / "cases").glob("*shadow.nc"))
    if len(nc_files) != 2:
        raise SystemExit(f"attesi due NetCDF shadow, trovati {len(nc_files)}")
    summaries = []
    for nc_path in nc_files:
        result = extract_case(nc_path, args.margin)
        write_case(result, output)
        summaries.append({key: value for key, value in result.items() if key != "nodal_loads"})
        print(
            f"[loads] n_nom={result['metadata']['nominal_load_factor']:.1f} "
            f"n_body={result['specific_acceleration_body_g_mean'][2]:.4f} "
            f"Fz={result['aerodynamic_force_body_lbf_mean'][2]:.3f} lbf"
        )
    (output / "summary.json").write_text(json.dumps(summaries, indent=2))
    (output / "run_all_zeno.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        "root=$(cd -- \"$(dirname -- \"${BASH_SOURCE[0]}\")\" && pwd)\n"
        "(cd \"$root/n1p0\" && ./run_supported_unloaded_zeno.sh)\n"
        "(cd \"$root/n1p0\" && ./run_zeno.sh)\n"
        "(cd \"$root/n1p6\" && ./run_zeno.sh)\n"
    )
    (output / "run_all_zeno.sh").chmod(0o755)
    print(f"[written] {output}")


if __name__ == "__main__":
    main()
