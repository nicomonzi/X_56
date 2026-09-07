#!/usr/bin/env python3
"""Analyse the completed 55 m/s coupled trim and create publication-ready plots."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file

WEIGHT_LBF = 419.4399625
FREEZE_TIME = 6.0


def slope(time: np.ndarray, values: np.ndarray) -> float:
    return float(np.polyfit(time, values, 1)[0])


def read_coupling(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {name: np.asarray([float(row[name]) for row in rows])
            for name in rows[0]}


def read_mbdyn(path: Path) -> dict[str, np.ndarray]:
    names = (
        "time", "run.iterations",
        "elem.joint.23.F", "elem.joint.23.M",
        "node.struct.990000.Phi",
        "elem.joint.1004.Phi", "elem.joint.2004.Phi",
        "elem.loadable.9101.output", "elem.loadable.9101.output_unsat",
        "elem.loadable.9101.error",
        "elem.loadable.9102.output", "elem.loadable.9102.output_unsat",
        "elem.loadable.9102.error", "elem.joint.5.a",
    )
    with netcdf_file(path, "r", mmap=False) as dataset:
        return {name: np.asarray(dataset.variables[name].data, dtype=float)
                for name in names}


def decorate(ax, ylabel: str) -> None:
    ax.axvspan(0.0, 0.4, color="0.92", label="rampa carichi")
    ax.axvline(0.4, color="0.45", lw=0.8, ls="--")
    ax.axvline(FREEZE_TIME, color="black", lw=1.0, ls=":")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    ax.set_xlim(0.0, 8.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path,
                        default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    work = args.case / "work" / "current"
    output = args.output or args.case / "results" / "analysis_55"
    output.mkdir(parents=True, exist_ok=True)

    aero = read_coupling(work / "coupling.csv")
    mb = read_mbdyn(work / "case.nc")
    ta = aero["time_s"]
    tm = mb["time"]
    af = ta >= FREEZE_TIME
    mf = tm >= FREEZE_TIME
    if not af.any() or not mf.any():
        raise RuntimeError("finestra congelata 6--8 s assente")

    pitch = np.rad2deg(mb["node.struct.990000.Phi"][:, 1])
    bfl = np.rad2deg(mb["elem.joint.1004.Phi"][:, 1])
    bfr = np.rad2deg(mb["elem.joint.2004.Phi"][:, 1])
    pitch_pid = np.rad2deg(mb["elem.loadable.9101.output"])
    flap_pid = np.rad2deg(mb["elem.loadable.9102.output"])
    joint_fz = mb["elem.joint.23.F"][:, 2]
    joint_my = mb["elem.joint.23.M"][:, 1]
    aero_fz = aero["total_fz_lbf"]
    aero_my = aero["global_my_cg_lbfin"]
    lift_excess = aero_fz - WEIGHT_LBF
    tip_l = aero["left_tip_dz_in"]
    tip_r = aero["right_tip_dz_in"]
    tip_s = aero["symmetric_tip_dz_in"]
    modal = mb["elem.joint.5.a"]

    saturated = np.flatnonzero(np.abs(bfl) >= 9.999)
    saturation_time = float(tm[saturated[0]]) if saturated.size else None
    modal_rows = {}
    for col, mode in enumerate(range(7, 19)):
        values = modal[mf, col]
        trend = np.polyval(np.polyfit(tm[mf], values, 1), tm[mf])
        modal_rows[str(mode)] = {
            "mean": float(np.mean(values)),
            "peak_to_peak": float(np.ptp(values)),
            "detrended_rms": float(np.sqrt(np.mean((values - trend) ** 2))),
        }

    summary = {
        "simulation_complete": bool(tm[-1] >= 7.996 - 1.e-6
                                    and ta[-1] >= 7.996 - 1.e-6),
        "last_time_s": float(min(tm[-1], ta[-1])),
        "trim_accepted": False,
        "reason": "BFL/BFR saturated at +10 deg with large nonzero Fz/My constraint reactions",
        "frozen_window_s": [6.0, float(min(tm[-1], ta[-1]))],
        "pitch_deg": {"mean": float(np.mean(pitch[mf])), "final": float(pitch[-1])},
        "bfl_deg": {"mean": float(np.mean(bfl[mf])), "final": float(bfl[-1])},
        "bfr_deg": {"mean": float(np.mean(bfr[mf])), "final": float(bfr[-1])},
        "body_flap_saturation_time_s": saturation_time,
        "aerodynamic_fz_lbf": {
            "mean": float(np.mean(aero_fz[af])), "final": float(aero_fz[-1]),
            "slope_lbf_s": slope(ta[af], aero_fz[af]),
        },
        "target_weight_lbf": WEIGHT_LBF,
        "aerodynamic_lift_excess_lbf": {
            "mean": float(np.mean(lift_excess[af])), "final": float(lift_excess[-1]),
        },
        "aerodynamic_my_lbfin": {
            "mean": float(np.mean(aero_my[af])), "final": float(aero_my[-1]),
            "slope_lbfin_s": slope(ta[af], aero_my[af]),
        },
        "trim_joint_fz_lbf": {
            "mean": float(np.mean(joint_fz[mf])), "final": float(joint_fz[-1]),
            "slope_lbf_s": slope(tm[mf], joint_fz[mf]),
        },
        "trim_joint_my_lbfin": {
            "mean": float(np.mean(joint_my[mf])), "final": float(joint_my[-1]),
            "slope_lbfin_s": slope(tm[mf], joint_my[mf]),
        },
        "tip_deflection_in": {
            "left_mean": float(np.mean(tip_l[af])),
            "right_mean": float(np.mean(tip_r[af])),
            "symmetric_mean": float(np.mean(tip_s[af])),
            "symmetric_final": float(tip_s[-1]),
            "symmetric_slope_in_s": slope(ta[af], tip_s[af]),
        },
        "solver_iterations": {
            "mean": float(np.mean(mb["run.iterations"])),
            "max": int(np.max(mb["run.iterations"])),
        },
        "modal_coordinates_frozen": modal_rows,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    plt.style.use("seaborn-v0_8-whitegrid")
    colors = plt.get_cmap("tab10").colors

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9), sharex=True,
                             layout="constrained")
    ax = axes[0, 0]
    ax.plot(tm, pitch, lw=1.8, label="pitch effettivo")
    ax.plot(tm, pitch_pid, lw=1.0, ls="--", label="uscita PID pitch")
    decorate(ax, "Pitch [deg]"); ax.legend(loc="best")
    ax = axes[0, 1]
    ax.plot(tm, bfl, lw=1.8, label="BFL")
    ax.plot(tm, bfr, lw=1.2, ls="--", label="BFR")
    ax.plot(tm, flap_pid, lw=1.0, ls=":", label="uscita PID flap")
    ax.axhline(10, color="tab:red", lw=0.8, ls="--", label="limite")
    decorate(ax, "Deflessione [deg]"); ax.legend(loc="best")
    ax = axes[1, 0]
    ax.plot(tm, joint_fz, lw=1.5, label="reazione Fz")
    ax.axhline(0, color="black", lw=0.8)
    decorate(ax, "Fz vincolo [lbf]"); ax.set_xlabel("Tempo [s]")
    ax = axes[1, 1]
    ax.plot(tm, joint_my, lw=1.5, color="tab:red", label="reazione My")
    ax.axhline(0, color="black", lw=0.8)
    decorate(ax, "My vincolo [lbf in]"); ax.set_xlabel("Tempo [s]")
    fig.suptitle("Trim online: comandi e residui del vincolo")
    fig.savefig(output / "trim_controls_and_residuals.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9.5), sharex=True,
                             layout="constrained")
    ax = axes[0, 0]
    ax.plot(ta, aero_fz, lw=1.6, label="Fz aerodinamica")
    ax.axhline(WEIGHT_LBF, color="black", lw=1.0, ls="--", label="peso target")
    decorate(ax, "Fz [lbf]"); ax.legend(loc="best")
    ax = axes[0, 1]
    ax.plot(ta, aero_my, lw=1.6, color="tab:red")
    ax.axhline(0, color="black", lw=0.8)
    decorate(ax, "My al CG [lbf in]")
    ax = axes[1, 0]
    ax.plot(ta, lift_excess, lw=1.5, label="Fz aero - peso")
    ax.plot(tm, joint_fz, lw=1.0, ls="--", label="reazione Fz")
    ax.axhline(0, color="black", lw=0.8)
    decorate(ax, "Squilibrio verticale [lbf]"); ax.set_xlabel("Tempo [s]"); ax.legend()
    ax = axes[1, 1]
    ax.plot(ta, aero_my, lw=1.5, label="My aero")
    ax.plot(tm, joint_my, lw=1.0, ls="--", label="reazione My")
    ax.axhline(0, color="black", lw=0.8)
    decorate(ax, "Momento [lbf in]"); ax.set_xlabel("Tempo [s]"); ax.legend()
    fig.suptitle("Bilancio aerodinamico e reazioni di trim")
    fig.savefig(output / "aerodynamic_balance.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                             layout="constrained")
    ax = axes[0]
    ax.plot(ta, tip_l, lw=1.2, label="tip sinistra")
    ax.plot(ta, tip_r, lw=1.2, label="tip destra")
    ax.plot(ta, tip_s, lw=1.8, color="black", label="media simmetrica")
    decorate(ax, "Deflessione tip [in]"); ax.legend(ncol=3)
    ax = axes[1]
    for col, mode in enumerate(range(7, 19)):
        ax.plot(tm, modal[:, col], lw=0.9, color=colors[col % 10],
                label=f"modo {mode}")
    decorate(ax, "Coordinata modale q"); ax.set_xlabel("Tempo [s]")
    ax.legend(ncol=6, fontsize=8)
    fig.suptitle("Risposta elastica")
    fig.savefig(output / "structural_response.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                             layout="constrained")
    ax = axes[0]
    ax.plot(tm, mb["elem.loadable.9101.error"], label="errore combinato pitch")
    ax.plot(tm, mb["elem.loadable.9102.error"], label="errore combinato flap")
    decorate(ax, "Ingresso PID"); ax.legend()
    ax = axes[1]
    ax.plot(tm, mb["run.iterations"], lw=0.9)
    decorate(ax, "Iterazioni MBDyn"); ax.set_xlabel("Tempo [s]")
    fig.suptitle("Attività dei controllori e convergenza numerica")
    fig.savefig(output / "controller_and_solver.png", dpi=180)
    plt.close(fig)

    report = f"""# Analisi trim accoppiato a 55 m/s

La simulazione è completa fino a {summary['last_time_s']:.3f} s ed è numericamente stabile,
ma **il trim non è stato raggiunto**.

- Pitch congelato: {summary['pitch_deg']['mean']:.4f} deg.
- BFL/BFR congelati: {summary['bfl_deg']['mean']:.4f} / {summary['bfr_deg']['mean']:.4f} deg;
  saturazione a +10 deg da t={saturation_time:.3f} s.
- Fz aerodinamica media 6--8 s: {summary['aerodynamic_fz_lbf']['mean']:.3f} lbf,
  contro {WEIGHT_LBF:.3f} lbf.
- Reazione verticale media residua: {summary['trim_joint_fz_lbf']['mean']:.3f} lbf.
- Momento aerodinamico medio al CG: {summary['aerodynamic_my_lbfin']['mean']:.3f} lbf in.
- Reazione di beccheggio media residua: {summary['trim_joint_my_lbfin']['mean']:.3f} lbf in.
- Deflessione simmetrica media delle estremità: {summary['tip_deflection_in']['symmetric_mean']:.4f} in.

Le pendenze nella finestra congelata sono piccole: il sistema si è assestato, ma
su un equilibrio imposto dai vincoli e non su una condizione di volo trimmata.
Il comando dei body flap procede nel verso che aumenta la deflessione positiva
fino alla saturazione senza annullare My: prima di una nuova run va corretto il
segno/il Jacobiano del canale BFL/BFR e va resa coerente la condizione iniziale.
Il seed usato ({pitch[100]:.3f} deg, {bfl[100]:.3f} deg circa) non è coerente
con il trim DUST precedentemente trovato vicino a 0.527 deg e -3.301 deg.

MBDyn converge con una media di {summary['solver_iterations']['mean']:.3f}
iterazioni e un massimo di {summary['solver_iterations']['max']}. Dopo il
congelamento non si osserva crescita delle coordinate modali; il modo 7 domina
la deformazione statica. Questo indica regolarità numerica del trim vincolato,
non costituisce una verifica di flutter o del SAS.
"""
    (output / "REPORT.md").write_text(report)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
