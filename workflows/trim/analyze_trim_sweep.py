#!/usr/bin/env python3
"""Analizza lo sweep del trim 2x2 e crea CSV e grafici riepilogativi."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mbdyn_trim_2x2_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.ticker import ScalarFormatter
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parent
DEFAULT_RESULTS = Path("/mnt/c/Users/Utente/Desktop/TRIM")
RAD_TO_DEG = 180.0 / math.pi
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167
IPS_DENSITY_TO_SI = 10686895.178201316
THESIS_COLORS = [
    "#8B0000",  # dark red
    "#00008B",  # dark blue
    "#4B0082",  # dark purple
    "#006400",  # dark green
    "#CC5500",  # dark orange
]
MULTI_TITLE_SIZE = 12
MULTI_LABEL_SIZE = 12
MULTI_LEGEND_SIZE = 12
MULTI_TICK_SIZE = 11
GRID_TITLE_SIZE = 13
GRID_LABEL_SIZE = 13
GRID_LEGEND_SIZE = 13
GRID_TICK_SIZE = 11

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.weight": "normal",
        "axes.titleweight": "normal",
        "axes.labelweight": "normal",
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "axes.prop_cycle": cycler(color=THESIS_COLORS),
        "lines.linewidth": 1.8,
    }
)


@dataclass
class Case:
    velocity: float
    density: float
    time: np.ndarray
    pitch: np.ndarray
    elevator: np.ndarray
    fz: np.ndarray
    my: np.ndarray
    pitch_pid: np.ndarray
    elevator_pid: np.ndarray
    directory: Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--fz-limit", type=float, default=5.0, help="soglia |Fz| [N]")
    parser.add_argument("--my-limit", type=float, default=0.5, help="soglia |My| [N m]")
    parser.add_argument(
        "--slope-limit",
        type=float,
        default=0.005,
        help="soglia deriva finale dei comandi [deg/s]",
    )
    return parser.parse_args()


def data(nc: Dataset, name: str) -> np.ndarray:
    if name not in nc.variables:
        raise KeyError(f"Variabile NetCDF assente: {name}")
    return np.ma.filled(nc.variables[name][:], np.nan).astype(float)


def source_constant(source: str, name: str) -> float:
    match = re.search(
        rf"set\s*:\s*(?:ifndef\s+)?const\s+real\s+{re.escape(name)}\s*=\s*"
        r"([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+\-]?\d+)?)\s*;",
        source,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError(f"Costante {name} non trovata nell'input del caso")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def load_case(case_directory: Path) -> Case:
    source = (case_directory / "case_input.mbd").read_text(encoding="utf-8")
    velocity = source_constant(source, "V_INF")
    density = source_constant(source, "RHO_AIR")
    with Dataset(case_directory / "result.nc") as nc:
        time = data(nc, "time")
        pitch = data(nc, "node.struct.990000.Phi")[:, 1] * RAD_TO_DEG
        # BFL is the left body flap; BFL/BFR receive the same symmetric trim command.
        elevator = data(nc, "elem.joint.1004.Phi")[:, 1] * RAD_TO_DEG
        fz = data(nc, "elem.joint.23.F")[:, 2] * LBF_TO_N
        my = data(nc, "elem.joint.23.M")[:, 1] * LBFIN_TO_NM
        pitch_pid = data(nc, "elem.loadable.9101.output").squeeze() * RAD_TO_DEG
        elevator_pid = data(nc, "elem.loadable.9102.output").squeeze() * RAD_TO_DEG
    count = min(map(len, (time, pitch, elevator, fz, my, pitch_pid, elevator_pid)))
    return Case(
        velocity,
        density,
        time[:count],
        pitch[:count],
        elevator[:count],
        fz[:count],
        my[:count],
        pitch_pid[:count],
        elevator_pid[:count],
        case_directory,
    )


def moving_average(values: np.ndarray, samples: int) -> np.ndarray:
    samples = max(1, int(samples))
    if samples == 1:
        return values.copy()
    left = samples // 2
    right = samples - 1 - left
    padded = np.pad(values, (left, right), mode="edge")
    return np.convolve(padded, np.ones(samples) / samples, mode="valid")


def convergence_time(case: Case, fz_limit: float, my_limit: float) -> float:
    dt = float(np.median(np.diff(case.time)))
    samples = max(3, int(round(1.0 / dt)))
    fz_slow = moving_average(case.fz, samples)
    my_slow = moving_average(case.my, samples)
    active = (case.time >= 2.0) & (case.time < 20.0)
    indices = np.flatnonzero(active)
    if not len(indices):
        return math.nan
    violation = (np.abs(fz_slow) > fz_limit) | (np.abs(my_slow) > my_limit)
    bad = indices[violation[indices]]
    if not len(bad):
        return float(case.time[indices[0]])
    candidate = int(bad[-1] + 1)
    if candidate >= len(case.time) or case.time[candidate] >= 20.0:
        return math.nan
    return float(case.time[candidate])


def slope(time: np.ndarray, values: np.ndarray, mask: np.ndarray) -> float:
    if np.count_nonzero(mask) < 2:
        return math.nan
    return float(np.polyfit(time[mask], values[mask], 1)[0])


def save_history(case: Case) -> None:
    destination = case.directory / "history.csv"
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "time_s",
                "pitch_deg",
                "elevator_deg",
                "Fz_N",
                "My_Nm",
                "pitch_pid_deg",
                "elevator_pid_deg",
            ]
        )
        writer.writerows(
            zip(
                case.time,
                case.pitch,
                case.elevator,
                case.fz,
                case.my,
                case.pitch_pid,
                case.elevator_pid,
            )
        )


def analyse(case: Case, args: argparse.Namespace) -> dict[str, object]:
    prefreeze = (case.time >= 18.0) & (case.time < 20.0)
    settled = case.time >= 22.0
    controlled = (case.time >= 2.0) & (case.time < 20.0)
    if np.count_nonzero(prefreeze) < 10 or np.count_nonzero(settled) < 10:
        raise RuntimeError(f"Risultato incompleto: tempo finale {case.time[-1]:.3f} s")

    pitch_slope = slope(case.time, case.pitch_pid, prefreeze)
    elevator_slope = slope(case.time, case.elevator_pid, prefreeze)
    fz_mean = float(np.nanmean(case.fz[settled]))
    my_mean = float(np.nanmean(case.my[settled]))
    saturated = bool(
        np.nanmax(np.abs(case.pitch_pid[controlled])) >= 14.99
        or np.nanmax(np.abs(case.elevator_pid[controlled])) >= 9.99
    )
    converged = bool(
        abs(fz_mean) <= args.fz_limit
        and abs(my_mean) <= args.my_limit
        and abs(pitch_slope) <= args.slope_limit
        and abs(elevator_slope) <= args.slope_limit
        and not saturated
    )
    rho_si = case.density * IPS_DENSITY_TO_SI
    return {
        "velocity_mps": case.velocity,
        "density_ips": case.density,
        "density_kg_m3": rho_si,
        "dynamic_pressure_Pa": 0.5 * rho_si * case.velocity**2,
        "pitch_trim_deg": float(np.nanmean(case.pitch[settled])),
        "elevator_trim_deg": float(np.nanmean(case.elevator[settled])),
        "Fz_mean_N": fz_mean,
        "Fz_std_N": float(np.nanstd(case.fz[settled])),
        "My_mean_Nm": my_mean,
        "My_std_Nm": float(np.nanstd(case.my[settled])),
        "pitch_slope_deg_s": pitch_slope,
        "elevator_slope_deg_s": elevator_slope,
        "convergence_time_s": convergence_time(case, args.fz_limit, args.my_limit),
        "pitch_pid_peak_deg": float(np.nanmax(np.abs(case.pitch_pid[controlled]))),
        "elevator_pid_peak_deg": float(np.nanmax(np.abs(case.elevator_pid[controlled]))),
        "saturated": saturated,
        "status": "converged" if converged else "check",
    }


def write_summary(rows: list[dict[str, object]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig: plt.Figure, output: Path, name: str) -> None:
    png_path = output / name
    pdf_path = output / f"{Path(name).stem}.pdf"
    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def format_axis(
    axis: plt.Axes,
    title: str,
    title_size: int,
    label_size: int,
    tick_size: int,
) -> None:
    axis.set_title(title, fontsize=title_size, fontweight="normal")
    axis.tick_params(
        axis="both",
        which="major",
        direction="in",
        labelsize=tick_size,
        length=5,
        width=0.8,
        top=True,
        right=True,
    )
    axis.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=3,
        width=0.6,
        top=True,
        right=True,
    )
    axis.minorticks_on()
    axis.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.35)
    axis.xaxis.label.set_size(label_size)
    axis.yaxis.label.set_size(label_size)
    axis.xaxis.label.set_fontweight("normal")
    axis.yaxis.label.set_fontweight("normal")
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))
    axis.yaxis.set_major_formatter(formatter)
    axis.yaxis.get_offset_text().set_fontsize(tick_size)
    axis.xaxis.get_offset_text().set_fontsize(tick_size)
    for label in (*axis.get_xticklabels(), *axis.get_yticklabels()):
        label.set_fontweight("normal")


def format_legend(axis: plt.Axes, size: int, **kwargs: object) -> None:
    legend = axis.legend(
        fontsize=size,
        frameon=True,
        framealpha=0.90,
        edgecolor="0.75",
        fancybox=False,
        **kwargs,
    )
    for text in legend.get_texts():
        text.set_fontweight("normal")


def plot_solution(rows: list[dict[str, object]], output: Path) -> None:
    v = np.asarray([row["velocity_mps"] for row in rows], dtype=float)
    pitch = np.asarray([row["pitch_trim_deg"] for row in rows], dtype=float)
    elevator = np.asarray([row["elevator_trim_deg"] for row in rows], dtype=float)
    plot_mask = ~np.isclose(v, 30.0)
    v_plot = v[plot_mask]
    pitch_plot = pitch[plot_mask]
    elevator_plot = elevator[plot_mask]
    fig, axes = plt.subplots(
        2, 1, figsize=(6.4, 7.2), sharex=True, constrained_layout=True
    )
    axes[0].plot(v_plot, pitch_plot, marker="o", linewidth=1.8, color="#00008B")
    axes[1].plot(
        v_plot, elevator_plot, marker="o", linewidth=1.8, color="#8B0000"
    )
    axes[0].set_ylabel(r"$\theta$ [deg]")
    axes[1].set_ylabel(r"$\delta_e$ [deg]")
    axes[1].set_xlabel(r"$V_\infty$ [m/s]")
    format_axis(
        axes[0],
        "Pitch Angle at Trim Condition",
        MULTI_TITLE_SIZE,
        MULTI_LABEL_SIZE,
        MULTI_TICK_SIZE,
    )
    format_axis(
        axes[1],
        "Elevator Deflection at Trim Condition",
        MULTI_TITLE_SIZE,
        MULTI_LABEL_SIZE,
        MULTI_TICK_SIZE,
    )
    save_figure(fig, output, "01_trim_solution.png")


def plot_quality(rows: list[dict[str, object]], output: Path) -> None:
    v = np.asarray([row["velocity_mps"] for row in rows], dtype=float)
    fz = np.asarray([row["Fz_mean_N"] for row in rows], dtype=float)
    fz_std = np.asarray([row["Fz_std_N"] for row in rows], dtype=float)
    my = np.asarray([row["My_mean_Nm"] for row in rows], dtype=float)
    my_std = np.asarray([row["My_std_Nm"] for row in rows], dtype=float)
    conv = np.asarray([row["convergence_time_s"] for row in rows], dtype=float)
    fig, axes = plt.subplots(
        3, 1, figsize=(6.4, 10.5), sharex=True, constrained_layout=True
    )
    axes[0].errorbar(v, fz, yerr=fz_std, marker="o", capsize=3, linewidth=1.8)
    axes[1].errorbar(v, my, yerr=my_std, marker="o", capsize=3, linewidth=1.8)
    axes[2].plot(v, conv, marker="o", linewidth=1.8)
    axes[0].axhline(0.0, color="0.2", linewidth=1.2, linestyle="-")
    axes[1].axhline(0.0, color="0.2", linewidth=1.2, linestyle="-")
    axes[0].set_ylabel(r"$F_z$ [N]")
    axes[1].set_ylabel(r"$M_y$ [N m]")
    axes[2].set_ylabel(r"$t_c$ [s]")
    axes[2].set_xlabel(r"$V_\infty$ [m/s]")
    titles = ["Vertical Force", "Pitching Moment", "Trim Convergence"]
    for axis, title in zip(axes, titles):
        format_axis(
            axis, title, MULTI_TITLE_SIZE, MULTI_LABEL_SIZE, MULTI_TICK_SIZE
        )
    save_figure(fig, output, "02_trim_quality.png")


def plot_histories(cases: list[Case], output: Path) -> None:
    targets = [35.0, 50.0, 70.0]
    selected = [min(cases, key=lambda case: abs(case.velocity - target)) for target in targets]
    fig, axes = plt.subplots(
        2, 2, figsize=(10.0, 8.0), sharex=True, constrained_layout=True
    )
    history_colors = ["#8B0000", "#00008B", "#56B4E9"]
    for case, color in zip(selected, history_colors):
        label = f"{case.velocity:g} m/s"
        plot_options = dict(
            label=label, linestyle="-", linewidth=1.8, color=color
        )
        axes[0, 0].plot(case.time, case.pitch, **plot_options)
        axes[0, 1].plot(case.time, case.elevator, **plot_options)
        axes[1, 0].plot(case.time, case.fz, **plot_options)
        axes[1, 1].plot(case.time, case.my, **plot_options)
    labels = [r"$\theta$ [deg]", r"$\delta_e$ [deg]", r"$F_z$ [N]", r"$M_y$ [N m]"]
    titles = ["Pitch angle", "Elevator Deflection", "Vertical Force", "Pitching Moment"]
    for axis, label, title in zip(axes.flat, labels, titles):
        axis.axvline(20.0, color="0.2", linestyle="-", linewidth=1.2)
        axis.set_ylabel(label)
        axis.set_xlabel(r"$t$ [s]")
        format_axis(axis, title, GRID_TITLE_SIZE, GRID_LABEL_SIZE, GRID_TICK_SIZE)
        format_legend(axis, GRID_LEGEND_SIZE, loc="best")
    save_figure(fig, output, "03_selected_histories.png")


def run() -> int:
    args = arguments()
    results = args.results.expanduser().resolve()
    output = (args.output.expanduser().resolve() if args.output else results / "analysis")
    output.mkdir(parents=True, exist_ok=True)

    cases: list[Case] = []
    failures: list[str] = []
    for case_directory in sorted(results.glob("V_*_mps")):
        if not (case_directory / "result.nc").is_file():
            failures.append(f"{case_directory.name}: result.nc assente")
            continue
        try:
            case = load_case(case_directory)
            save_history(case)
            cases.append(case)
        except Exception as error:  # mostra tutti i casi problematici nel report
            failures.append(f"{case_directory.name}: {error}")
    cases.sort(key=lambda case: case.velocity)
    if not cases:
        raise RuntimeError(f"Nessun caso NetCDF leggibile sotto {results}")

    rows = [analyse(case, args) for case in cases]
    write_summary(rows, output / "trim_summary.csv")
    plot_solution(rows, output)
    plot_quality(rows, output)
    plot_histories(cases, output)

    with (output / "analysis_report.txt").open("w", encoding="utf-8") as stream:
        stream.write("Sweep trim longitudinale accoppiato 2x2\n")
        stream.write(f"Casi analizzati: {len(rows)}\n")
        stream.write(f"Casi convergenti: {sum(row['status'] == 'converged' for row in rows)}\n")
        if failures:
            stream.write("\nCasi mancanti o illeggibili:\n")
            stream.write("\n".join(failures) + "\n")

    print("V [m/s]   pitch [deg]   elev [deg]    Fz [N]    My [Nm]   stato")
    for row in rows:
        print(
            f"{row['velocity_mps']:7.1f} {row['pitch_trim_deg']:13.5f} "
            f"{row['elevator_trim_deg']:12.5f} {row['Fz_mean_N']:9.3f} "
            f"{row['My_mean_Nm']:9.4f}   {row['status']}"
        )
    print(f"\nRisultati analisi: {output}")
    if failures:
        print(f"Attenzione: {len(failures)} casi mancanti o illeggibili")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
