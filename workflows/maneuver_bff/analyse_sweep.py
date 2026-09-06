#!/usr/bin/env python3
"""Analyze paired shadow/excited results with one method for all campaigns."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/manouver_stifness_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from campaign import CONFIG, MANEUVER, MANEUVER_DIR

if str(MANEUVER_DIR) not in sys.path:
    sys.path.insert(0, str(MANEUVER_DIR))
from analyse_time_domain_pairs import (  # noqa: E402
    constant,
    paired_body_flap_leakage,
    paired_tip_delta,
    trajectory_metrics,
)
from compare_paired_response import growth_metrics, paired_delta  # noqa: E402

IN_TO_M = 0.0254
META_RE = re.compile(r"(?m)^# STIFFNESS_STUDY_METADATA (\{.*\})$")


def metadata(path: Path) -> dict:
    match = META_RE.search(path.read_text())
    if not match:
        raise ValueError(f"metadata stiffness assente in {path}")
    return json.loads(match.group(1))


def is_complete(nc: Path, mbd: Path) -> bool:
    with Dataset(nc) as data:
        time = np.asarray(data["time"][:]).squeeze()
    required = constant(mbd, "SAS_OFF_START") + constant(mbd, "SAS_OFF_DURATION")
    return time.ndim == 1 and len(time) >= 50 and float(time[-1]) >= required - 0.01


def discover(case_dir: Path) -> dict[tuple, dict[bool, tuple[Path, Path, dict]]]:
    groups: dict[tuple, dict[bool, tuple[Path, Path, dict]]] = {}
    for nc in sorted(case_dir.glob("*.nc")):
        mbd = nc.with_suffix(".mbd")
        if not mbd.is_file():
            continue
        try:
            meta = metadata(mbd)
            if not is_complete(nc, mbd):
                print(f"[skip incomplete] {nc.name}")
                continue
        except Exception as error:
            print(f"[skip unreadable] {nc.name}: {error}")
            continue
        key = (
            str(meta["campaign"]), float(meta["velocity_mps"]),
            float(meta["nominal_load_factor"]), float(meta["time_step_s"]),
            float(meta["mode7_frequency_scale"]),
        )
        excited = bool(meta["excited"])
        if excited in groups.setdefault(key, {}):
            raise RuntimeError(f"duplicato per {key}, excited={excited}")
        groups[key][excited] = (nc, mbd, meta)
    return groups


def tas_metrics(nc: Path, mbd: Path) -> dict[str, float]:
    with Dataset(nc) as data:
        time = np.asarray(data["time"][:]).squeeze()
        velocity = np.asarray(data["node.struct.990000.XP"][:]) * IN_TO_M
    start = constant(mbd, "SAS_OFF_START")
    end = start + constant(mbd, "SAS_OFF_DURATION")
    use = (time >= start) & (time <= end)
    # XP is the aircraft perturbation velocity, while V_INF is the ambient
    # air velocity used by the aerodynamic elements. Their relative velocity
    # gives TAS; norm(XP) alone is not airspeed.
    air_velocity = np.array([constant(mbd, "V_INF"), 0.0, 0.0])
    tas = np.linalg.norm(air_velocity - velocity[use], axis=1)
    slope = np.polyfit(time[use] - time[use][0], tas, 1)[0]
    return {
        "tas_mean_sas_off_mps": float(np.mean(tas)),
        "tas_std_sas_off_mps": float(np.std(tas)),
        "tas_slope_sas_off_mps2": float(slope),
        "tas_min_sas_off_mps": float(np.min(tas)),
        "tas_max_sas_off_mps": float(np.max(tas)),
    }


def analyze_pairs(campaign_dir: Path) -> tuple[list[dict], list[dict]]:
    groups = discover(campaign_dir / "cases")
    acceptance = CONFIG["acceptance"]
    rows: list[dict] = []
    missing: list[dict] = []
    for key, pair in sorted(groups.items()):
        campaign, velocity, nominal_n, dt, scale = key
        if set(pair) != {False, True}:
            missing.append({
                "campaign": campaign, "velocity_mps": velocity,
                "nominal_load_factor": nominal_n, "time_step_s": dt,
                "mode7_frequency_scale": scale,
                "missing": "excited" if False in pair else "shadow",
            })
            continue
        shadow_nc, shadow_mbd, shadow_meta = pair[False]
        excited_nc, excited_mbd, _ = pair[True]
        time, delta_q7 = paired_delta(shadow_nc, excited_nc)
        _, delta_tip = paired_tip_delta(shadow_nc, excited_nc)
        frequency = float(MANEUVER.baseline_module().dlm_rom_values(velocity)[3])
        release = constant(excited_mbd, "SAS_OFF_START")
        duration = constant(excited_mbd, "SAS_OFF_DURATION")
        rap_duration = 0.742 / frequency
        fit_start = release + 0.05 + rap_duration + 0.05
        fit_end = release + duration - 0.05
        mode_growth, _, _ = growth_metrics(
            time, delta_q7, fit_start, fit_end, frequency
        )
        tip_growth, _, _ = growth_metrics(
            time, delta_tip, fit_start, fit_end, frequency
        )
        trajectory = trajectory_metrics(shadow_nc, shadow_mbd)
        tas = tas_metrics(shadow_nc, shadow_mbd)
        leakage = paired_body_flap_leakage(
            shadow_nc, excited_nc, fit_start, fit_end, frequency
        )
        n_error = trajectory["achieved_n_mean_sas_off"] - nominal_n
        q_command = float(shadow_meta["pitch_rate_command_deg_s"])
        q_error = trajectory["q_mean_sas_off_deg_s"] - q_command
        q_tolerance = max(
            float(acceptance["maximum_pitch_rate_absolute_error_deg_s"]),
            float(acceptance["maximum_pitch_rate_relative_error"]) * abs(q_command),
        )
        stationarity_valid = bool(
            abs(n_error) <= float(acceptance["maximum_abs_nominal_n_error"])
            and trajectory["achieved_n_std_sas_off"]
            <= float(acceptance["maximum_n_standard_deviation"])
            and abs(trajectory["achieved_n_slope_sas_off_per_s"])
            <= float(acceptance["maximum_abs_n_slope_per_s"])
            and tas["tas_std_sas_off_mps"]
            <= float(acceptance["maximum_tas_standard_deviation_mps"])
            and abs(tas["tas_mean_sas_off_mps"] - velocity)
            <= float(acceptance["maximum_tas_mean_error_mps"])
            and abs(q_error) <= q_tolerance
            and trajectory["q_std_sas_off_deg_s"]
            <= float(acceptance["maximum_pitch_rate_standard_deviation_deg_s"])
        )
        observable_consistency = abs(
            mode_growth["sigma_per_s"] - tip_growth["sigma_per_s"]
        )
        identification_valid = bool(
            stationarity_valid
            and observable_consistency
            <= float(acceptance["maximum_mode_tip_sigma_difference_per_s"])
            and leakage <= float(acceptance["maximum_surface_bff_leakage_deg"])
        )
        sigma = mode_growth["sigma_per_s"]
        resolution = float(CONFIG["physics"]["sigma_resolution_per_s"])
        rows.append({
            "campaign": campaign,
            "velocity_mps": velocity,
            "nominal_load_factor": nominal_n,
            "time_step_s": dt,
            "mode7_frequency_scale": scale,
            "mode7_frequency_shift_percent": 100.0 * (scale - 1.0),
            "fit_start_s": fit_start,
            "fit_end_s": fit_end,
            **trajectory,
            **tas,
            "nominal_n_error": n_error,
            "pitch_rate_command_deg_s": q_command,
            "pitch_rate_error_deg_s": q_error,
            "pitch_rate_tolerance_deg_s": q_tolerance,
            "mode7_sigma_per_s": sigma,
            "tip_sigma_per_s": tip_growth["sigma_per_s"],
            "mode_tip_abs_sigma_difference_per_s": observable_consistency,
            "surface_bff_leakage_deg": leakage,
            "stationarity_valid": stationarity_valid,
            "identification_valid": identification_valid,
            "verdict": (
                "invalid" if not identification_valid
                else "unstable" if sigma > resolution
                else "stable" if sigma < -resolution
                else "near_onset"
            ),
        })
    return rows, missing


def onset_summary(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        if row["identification_valid"]:
            grouped[(
                row["campaign"], row["nominal_load_factor"], row["time_step_s"],
                row["mode7_frequency_scale"],
            )].append(row)
    result: list[dict] = []
    for key, points in sorted(grouped.items()):
        points = sorted(points, key=lambda item: item["velocity_mps"])
        bracket = None
        estimate = None
        for left, right in zip(points, points[1:]):
            sl, sr = left["mode7_sigma_per_s"], right["mode7_sigma_per_s"]
            if sl == 0.0 or sl * sr <= 0.0:
                vl, vr = left["velocity_mps"], right["velocity_mps"]
                bracket = (vl, vr)
                estimate = vl if sl == 0.0 else vl - sl * (vr - vl) / (sr - sl)
                break
        result.append({
            "campaign": key[0], "nominal_load_factor": key[1],
            "time_step_s": key[2], "mode7_frequency_scale": key[3],
            "valid_point_count": len(points),
            "status": "crossing_found" if bracket else "no_bracket",
            "lower_velocity_mps": None if bracket is None else bracket[0],
            "upper_velocity_mps": None if bracket is None else bracket[1],
            "onset_velocity_linear_mps": estimate,
            "mean_achieved_n": float(np.mean([
                point["achieved_n_mean_sas_off"] for point in points
            ])),
        })
    return result


def timestep_summary(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["velocity_mps"], row["nominal_load_factor"],
                 row["mode7_frequency_scale"])].append(row)
    result = []
    for key, points in sorted(grouped.items()):
        by_dt = {point["time_step_s"]: point for point in points}
        if 0.01 not in by_dt or 0.005 not in by_dt:
            continue
        coarse, fine = by_dt[0.01], by_dt[0.005]
        ds = fine["mode7_sigma_per_s"] - coarse["mode7_sigma_per_s"]
        dn = fine["achieved_n_mean_sas_off"] - coarse["achieved_n_mean_sas_off"]
        result.append({
            "velocity_mps": key[0], "nominal_load_factor": key[1],
            "mode7_frequency_scale": key[2],
            "sigma_dt_0p01_per_s": coarse["mode7_sigma_per_s"],
            "sigma_dt_0p005_per_s": fine["mode7_sigma_per_s"],
            "delta_sigma_fine_minus_coarse_per_s": ds,
            "delta_n_fine_minus_coarse": dn,
            "converged": bool(
                abs(ds) <= float(CONFIG["acceptance"]["maximum_timestep_sigma_difference_per_s"])
                and abs(dn) <= float(CONFIG["acceptance"]["maximum_timestep_n_difference"])
            ),
        })
    return result


def stiffness_summary(rows: list[dict]) -> tuple[list[dict], dict]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        if row["identification_valid"]:
            grouped[(row["velocity_mps"], row["nominal_load_factor"],
                     row["time_step_s"])].append(row)
    result = []
    predicted_spans = []
    for key, points in sorted(grouped.items()):
        if len({point["mode7_frequency_scale"] for point in points}) < 3:
            continue
        points = sorted(points, key=lambda point: point["mode7_frequency_scale"])
        x = np.asarray([point["mode7_frequency_scale"] - 1.0 for point in points])
        y = np.asarray([point["mode7_sigma_per_s"] for point in points])
        slope, intercept = np.polyfit(x, y, 1)
        residual = y - (slope * x + intercept)
        predicted_span = abs(float(slope)) * 0.02  # from -1% to +1%
        predicted_spans.append(predicted_span)
        result.append({
            "velocity_mps": key[0], "nominal_load_factor": key[1],
            "time_step_s": key[2], "valid_scale_count": len(points),
            "d_sigma_per_unit_frequency_scale": float(slope),
            "d_sigma_per_one_percent_frequency_shift_per_s": float(slope) * 0.01,
            "predicted_sigma_span_for_plus_minus_one_percent_per_s": predicted_span,
            "linear_fit_rms_per_s": float(np.sqrt(np.mean(residual * residual))),
            "measured_sigma_min_per_s": float(np.min(y)),
            "measured_sigma_max_per_s": float(np.max(y)),
        })
    resolution = float(CONFIG["physics"]["sigma_resolution_per_s"])
    if not predicted_spans:
        decision = "insufficient_results"
    elif max(predicted_spans) > resolution:
        decision = "physical_stress_stiffening_rom_recommended"
    else:
        decision = "prestress_not_resolved_skip_for_primary_conclusion"
    return result, {
        "decision": decision,
        "criterion": (
            "maximum predicted sigma span over a +/-1% mode-7 frequency shift "
            f"compared with {resolution:g} 1/s identification resolution"
        ),
        "maximum_predicted_span_per_s": None if not predicted_spans else max(predicted_spans),
        "warning": (
            "This is a parametric decision gate, not a prediction of the actual "
            "prestress-induced frequency shift."
        ),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_summary(rows: list[dict], output: Path) -> None:
    valid = [row for row in rows if row["identification_valid"]]
    if not valid:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    primary = [
        row for row in valid
        if row["campaign"] == "primary"
        and abs(row["mode7_frequency_scale"] - 1.0) < 1e-12
    ]
    for load in sorted({row["nominal_load_factor"] for row in primary}):
        group = sorted(
            (row for row in primary if row["nominal_load_factor"] == load),
            key=lambda row: row["velocity_mps"],
        )
        axes[0].plot(
            [row["velocity_mps"] for row in group],
            [row["mode7_sigma_per_s"] for row in group], "o-", label=f"n_nom={load:g}"
        )
    sensitivity_keys = {
        (row["velocity_mps"], row["nominal_load_factor"], row["time_step_s"])
        for row in valid
        if abs(row["mode7_frequency_scale"] - 1.0) > 1e-12
    }
    sensitivity = [
        row for row in valid
        if (row["velocity_mps"], row["nominal_load_factor"], row["time_step_s"])
        in sensitivity_keys
    ]
    for load in sorted({row["nominal_load_factor"] for row in sensitivity}):
        group = sorted(
            (row for row in sensitivity if row["nominal_load_factor"] == load),
            key=lambda row: row["mode7_frequency_scale"],
        )
        if len(group) > 1:
            axes[1].plot(
                [row["mode7_frequency_shift_percent"] for row in group],
                [row["mode7_sigma_per_s"] for row in group], "o-", label=f"n_nom={load:g}"
            )
    axes[0].set(xlabel="TAS [m/s]", ylabel="sigma [1/s]", title="Paired onset map")
    axes[1].set(xlabel="mode-7 frequency shift [%]", ylabel="sigma [1/s]",
                title="Parametric stiffness screen")
    for axis in axes:
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.grid(True, alpha=0.3)
        if axis.lines:
            axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "summary.png", dpi=190)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign_directory", type=Path)
    parser.add_argument(
        "--reference-directory", type=Path, action="append", default=[],
        help="completed campaign(s) supplying shared reference points",
    )
    args = parser.parse_args()
    campaign_dir = args.campaign_directory.expanduser().resolve()
    output = campaign_dir / "analysis"
    output.mkdir(parents=True, exist_ok=True)
    rows, missing = analyze_pairs(campaign_dir)
    for reference in args.reference_directory:
        extra_rows, extra_missing = analyze_pairs(reference.expanduser().resolve())
        rows.extend(extra_rows)
        missing.extend(extra_missing)
    if not rows and not missing:
        raise SystemExit("nessun risultato NetCDF completo trovato")
    onsets = onset_summary(rows)
    timestep = timestep_summary(rows)
    sensitivity, decision = stiffness_summary(rows)
    write_csv(output / "paired_results.csv", rows)
    write_csv(output / "onset.csv", onsets)
    write_csv(output / "timestep_convergence.csv", timestep)
    write_csv(output / "stiffness_sensitivity.csv", sensitivity)
    summary = {
        "target_campaign_directory": str(campaign_dir),
        "reference_directories": [
            str(path.expanduser().resolve()) for path in args.reference_directory
        ],
        "complete_pairs": len(rows),
        "valid_pairs": sum(bool(row["identification_valid"]) for row in rows),
        "missing_pairs": missing,
        "onsets": onsets,
        "timestep_convergence": timestep,
        "stiffness_sensitivity": sensitivity,
        "prestress_decision_gate": decision,
        "sigma_resolution_per_s": CONFIG["physics"]["sigma_resolution_per_s"],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    plot_summary(rows, output)
    print(json.dumps({
        "analysis": str(output), "complete_pairs": len(rows),
        "valid_pairs": summary["valid_pairs"],
        "prestress_decision": decision["decision"],
    }, indent=2))


if __name__ == "__main__":
    main()
