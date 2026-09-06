#!/usr/bin/env python3
"""Build a multimode root-locus diagram from the MBDyn free responses."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mbdyn_root_locus_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


DEFAULT_RESULTS = Path(
    "/mnt/c/Users/Utente/Desktop/RESULTS_BBF_DT002_MODES_7_23_DENSE"
)
MODAL_JOINT = 5
BASE_NODE = 990000
LEFT_TIP = 990020
RIGHT_TIP = 991020
INCH_TO_M = 0.0254
FREE_RESPONSE_DELAY_S = 0.50
MAX_FREQUENCY_HZ = 24.5

BRANCH_NAMES = (
    "Short-period",
    "Symmetric bending",
    "Symmetric torsion",
    "Antisymmetric bending",
    "Antisymmetric torsion",
)
# Finestre ampie centrate sui cinque rami osservati alla velocita' piu' bassa.
# Evitano che, al variare di V, due modi vicini si scambino semplicemente
# perche' cambia la loro ampiezza relativa nella risposta.
BRANCH_FREQUENCY_WINDOWS_HZ = (
    (0.20, 2.00),
    (1.80, 5.00),
    (5.00, 8.50),
    (8.20, 9.70),
    (9.40, 11.50),
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Identify the five principal rigid/elastic branches from every "
            "complete free response and draw a complex-plane root locus."
        )
    )
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Default: RESULTS/ROOT_LOCUS_MULTIMODE.",
    )
    return parser.parse_args()


def nc_array(dataset: Dataset, name: str) -> np.ndarray:
    if name not in dataset.variables:
        raise KeyError(f"Missing NetCDF variable: {name}")
    return np.ma.filled(dataset.variables[name][:], np.nan).astype(float)


def log_constant(path: Path, name: str, fallback: float) -> float:
    if not path.is_file():
        return fallback
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        rf"(?m)^\s*const\s+real\s+{re.escape(name)}\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s*$",
        text,
    )
    return fallback if match is None else float(match.group(1).replace("D", "E"))


def identify_poles(
    time: np.ndarray,
    generalized_signals: np.ndarray,
) -> list[dict[str, float | int | np.ndarray]]:
    """Identify continuous-time poles with a scaled full-state DMD operator."""
    derivatives = np.gradient(generalized_signals, time, axis=0)
    state = np.column_stack((generalized_signals, derivatives))
    state -= np.mean(state, axis=0, keepdims=True)
    scale = np.std(state, axis=0)
    # Un segnale non eccitato non deve generare un ramo basato sul rumore
    # numerico. La soglia e' applicata separatamente alle cinque misure.
    signal_scale = np.std(generalized_signals, axis=0)
    observable = signal_scale > max(
        1.0e-11, 1.0e-8 * float(np.max(signal_scale))
    )
    for branch_index, is_observable in enumerate(observable):
        if not is_observable:
            scale[branch_index] = np.nan
            scale[branch_index + len(BRANCH_NAMES)] = np.nan
    scale[scale <= 1.0e-14] = np.nan
    valid_channels = np.isfinite(scale)
    scaled = (state[:, valid_channels] / scale[valid_channels]).T

    x0 = scaled[:, :-1]
    x1 = scaled[:, 1:]
    u, singular_values, vh = np.linalg.svd(x0, full_matrices=False)
    if not len(singular_values) or singular_values[0] <= 0.0:
        return []
    rank = int(np.count_nonzero(singular_values >= singular_values[0] * 1.0e-8))
    rank = max(2, min(rank, 2 * len(BRANCH_NAMES)))
    u = u[:, :rank]
    singular_values = singular_values[:rank]
    vh = vh[:rank]

    # Operatore dinamico ridotto. Non viene eseguito alcun fit esponenziale
    # dell'inviluppo e non viene usato un criterio R^2.
    reduced_operator = (
        u.T @ x1 @ vh.T @ np.diag(1.0 / singular_values)
    )
    discrete_poles, reduced_shapes = np.linalg.eig(reduced_operator)
    shapes_scaled = u @ reduced_shapes
    dt = float(np.median(np.diff(time)))
    continuous_poles = np.log(discrete_poles.astype(complex)) / dt

    original_indices = np.flatnonzero(valid_channels)
    position_indices: list[int | None] = []
    velocity_indices: list[int | None] = []
    for branch_index in range(len(BRANCH_NAMES)):
        position_match = np.flatnonzero(original_indices == branch_index)
        velocity_match = np.flatnonzero(
            original_indices == branch_index + len(BRANCH_NAMES)
        )
        position_indices.append(
            int(position_match[0]) if len(position_match) else None
        )
        velocity_indices.append(
            int(velocity_match[0]) if len(velocity_match) else None
        )

    poles: list[dict[str, float | int | np.ndarray]] = []
    for pole_index, pole in enumerate(continuous_poles):
        frequency = float(pole.imag / (2.0 * math.pi))
        real_hz = float(pole.real / (2.0 * math.pi))
        if not (0.10 < frequency <= MAX_FREQUENCY_HZ):
            continue
        if not (-6.0 <= real_hz <= 3.0):
            continue

        shape = np.abs(shapes_scaled[:, pole_index])
        branch_participation = np.zeros(len(BRANCH_NAMES))
        for branch_index, (position_index, velocity_index) in enumerate(
            zip(position_indices, velocity_indices)
        ):
            if position_index is not None:
                branch_participation[branch_index] += shape[position_index] ** 2
            if velocity_index is not None:
                branch_participation[branch_index] += shape[velocity_index] ** 2
        total = float(np.sum(branch_participation))
        if total <= 0.0:
            continue
        poles.append(
            {
                "real_hz": real_hz,
                "frequency_hz": frequency,
                "branch_participation": branch_participation / total,
            }
        )
    return poles


def analyze_case(path: Path) -> dict[str, object]:
    log = path.with_suffix(".log")
    velocity = log_constant(log, "V_INF", math.nan)
    trim_end = log_constant(log, "TRIM_END", 15.0)
    excitation_off = log_constant(log, "BFF_WINDOW_START", 16.75)
    observation_end = log_constant(log, "BFF_WINDOW_END", 31.75)

    with Dataset(path) as dataset:
        time = nc_array(dataset, "time")
        base_position = nc_array(dataset, f"node.struct.{BASE_NODE}.X")
        base_attitude = nc_array(dataset, f"node.struct.{BASE_NODE}.Phi")
        left_position = nc_array(dataset, f"node.struct.{LEFT_TIP}.X")
        right_position = nc_array(dataset, f"node.struct.{RIGHT_TIP}.X")
        left_attitude = nc_array(dataset, f"node.struct.{LEFT_TIP}.Phi")
        right_attitude = nc_array(dataset, f"node.struct.{RIGHT_TIP}.Phi")
        modal_q = nc_array(dataset, f"elem.joint.{MODAL_JOINT}.a")

    if modal_q.ndim == 1:
        modal_q = modal_q[:, np.newaxis]
    complete = bool(len(time) and time[-1] >= observation_end)
    if modal_q.shape[1] != 17 or not complete:
        return {"velocity": velocity, "accepted": False, "poles": []}

    trim = (time >= trim_end - 1.0) & (time < trim_end)
    if np.count_nonzero(trim) < 10:
        return {"velocity": velocity, "accepted": False, "poles": []}

    # Deformazioni verticali delle tip depurate da heave e pitch rigidi.
    pitch = base_attitude[:, 1]
    pitch_relative = pitch - float(np.mean(pitch[trim]))
    z_cg_relative = (
        base_position[:, 2] - float(np.mean(base_position[trim, 2]))
    )
    x_cg_trim = float(np.mean(base_position[trim, 0]))
    x_left = float(np.mean(left_position[trim, 0])) - x_cg_trim
    x_right = float(np.mean(right_position[trim, 0])) - x_cg_trim
    rigid_left_z = z_cg_relative - x_left * pitch_relative
    rigid_right_z = z_cg_relative - x_right * pitch_relative
    elastic_left_z = (
        left_position[:, 2]
        - float(np.mean(left_position[trim, 2]))
        - rigid_left_z
    ) * INCH_TO_M
    elastic_right_z = (
        right_position[:, 2]
        - float(np.mean(right_position[trim, 2]))
        - rigid_right_z
    ) * INCH_TO_M

    # Le rotazioni elastiche sono relative al corpo rigido; somma e differenza
    # separano in modo esplicito i contributi simmetrici e antisimmetrici.
    left_elastic_pitch = left_attitude[:, 1] - pitch
    right_elastic_pitch = right_attitude[:, 1] - pitch
    generalized_signals = np.column_stack(
        (
            pitch_relative,
            0.5 * (elastic_left_z + elastic_right_z),
            0.5 * (left_elastic_pitch + right_elastic_pitch),
            0.5 * (elastic_left_z - elastic_right_z),
            0.5 * (left_elastic_pitch - right_elastic_pitch),
        )
    )

    free = (
        (time >= excitation_off + FREE_RESPONSE_DELAY_S)
        & (time <= observation_end)
    )
    finite = (
        np.isfinite(time)
        & np.all(np.isfinite(generalized_signals), axis=1)
    )
    selected = free & finite
    if np.count_nonzero(selected) < 200:
        return {"velocity": velocity, "accepted": False, "poles": []}

    poles = identify_poles(time[selected], generalized_signals[selected])
    return {"velocity": velocity, "accepted": bool(poles), "poles": poles}


def assign_branches(
    case: dict[str, object],
    previous_frequency: dict[str, float],
) -> list[dict[str, float | int | str]]:
    """Assign at most one identified pole to each physical output branch."""
    velocity = float(case["velocity"])
    poles = list(case["poles"])
    rows: list[dict[str, float | int | str]] = []

    for branch_index, branch_name in enumerate(BRANCH_NAMES):
        lower_frequency, upper_frequency = (
            BRANCH_FREQUENCY_WINDOWS_HZ[branch_index]
        )
        candidates = [
            pole
            for pole in poles
            if lower_frequency
            <= float(pole["frequency_hz"])
            <= upper_frequency
        ]
        if not candidates:
            continue
        if branch_name in previous_frequency:
            window_width = upper_frequency - lower_frequency
            pole = min(
                candidates,
                key=lambda item: (
                    abs(
                        float(item["frequency_hz"])
                        - previous_frequency[branch_name]
                    )
                    / window_width
                    - 0.35
                    * float(item["branch_participation"][branch_index])
                ),
            )
        else:
            pole = max(
                candidates,
                key=lambda item: float(
                    item["branch_participation"][branch_index]
                ),
            )
        participation = float(pole["branch_participation"][branch_index])
        if participation < 0.025:
            continue
        previous_frequency[branch_name] = float(pole["frequency_hz"])
        rows.append(
            {
                "velocity_mps": velocity,
                "branch": branch_name,
                "real_part_over_2pi_hz": float(pole["real_hz"]),
                "imaginary_part_over_2pi_hz": float(pole["frequency_hz"]),
                "participation": participation,
            }
        )
    return rows


def plot_root_locus(rows: list[dict[str, float | int | str]], path: Path) -> None:
    figure, axis = plt.subplots(figsize=(13.5, 8.0))
    velocities = np.array([float(row["velocity_mps"]) for row in rows])
    norm = plt.Normalize(float(np.min(velocities)), float(np.max(velocities)))
    colour_map = plt.get_cmap("viridis")

    available = {str(row["branch"]) for row in rows}
    branches = [name for name in BRANCH_NAMES if name in available]
    for branch_index, branch in enumerate(branches):
        branch_rows = sorted(
            (row for row in rows if row["branch"] == branch),
            key=lambda row: float(row["velocity_mps"]),
        )
        x = np.array(
            [float(row["real_part_over_2pi_hz"]) for row in branch_rows]
        )
        y = np.array(
            [float(row["imaginary_part_over_2pi_hz"]) for row in branch_rows]
        )
        v = np.array([float(row["velocity_mps"]) for row in branch_rows])
        line_colour = plt.cm.tab20(branch_index % 20)
        axis.plot(x, y, "-", color=line_colour, linewidth=0.9, alpha=0.65)
        axis.scatter(
            x,
            y,
            c=v,
            cmap=colour_map,
            norm=norm,
            s=29,
            edgecolors=line_colour,
            linewidths=0.55,
            zorder=3,
        )
        if len(x):
            label_index = int(np.argmin(x))
            axis.annotate(
                branch,
                xy=(x[label_index], y[label_index]),
                xytext=(-6, 5),
                textcoords="offset points",
                fontsize=7.5,
                color=line_colour,
                ha="right",
            )

    axis.axvline(0.0, color="black", linewidth=1.1)
    axis.text(
        0.99,
        0.985,
        "Unstable half-plane",
        transform=axis.transAxes,
        ha="right",
        va="top",
        color="tab:red",
        fontsize=9,
    )
    scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=colour_map)
    colour_bar = figure.colorbar(scalar_map, ax=axis, pad=0.015)
    colour_bar.set_label("Equivalent airspeed [m/s]")
    axis.set_xlabel(r"Real part, $\mathrm{Re}(\lambda)/(2\pi)$ [Hz]")
    axis.set_ylabel(r"Imaginary part, $\mathrm{Im}(\lambda)/(2\pi)$ [Hz]")
    axis.set_title(
        "Root loci of the principal rigid and elastic modes"
    )
    # Inquadra soltanto i cinque rami richiesti, evitando lo spazio vuoto fino
    # alla frequenza di Nyquist dei modi FEM che non vengono rappresentati.
    maximum_plotted_frequency = max(
        float(row["imaginary_part_over_2pi_hz"]) for row in rows
    )
    axis.set_ylim(0.0, math.ceil(maximum_plotted_frequency + 0.8))
    axis.grid(True, alpha=0.25)
    figure.tight_layout()
    figure.savefig(path, dpi=190)
    plt.close(figure)


def main() -> None:
    args = arguments()
    results = args.results.expanduser().resolve()
    output = (args.output or results / "ROOT_LOCUS_MULTIMODE").resolve()
    output.mkdir(parents=True, exist_ok=True)

    paths = sorted(results.glob("V_*_mps/case.nc"))
    if not paths:
        raise FileNotFoundError(f"No case.nc files found below {results}")
    cases = [analyze_case(path) for path in paths]
    cases = [case for case in cases if bool(case["accepted"])]
    cases.sort(key=lambda case: float(case["velocity"]))
    previous_frequency: dict[str, float] = {}
    rows = []
    for case in cases:
        rows.extend(assign_branches(case, previous_frequency))
    if not rows:
        raise RuntimeError("No observable oscillatory poles were identified")

    csv_path = output / "root_locus_poles.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure_path = output / "root_locus_multimode.png"
    plot_root_locus(rows, figure_path)

    notes_path = output / "root_locus_notes.txt"
    with notes_path.open("w", encoding="utf-8") as stream:
        stream.write("RESPONSE-IDENTIFIED MULTIMODE ROOT LOCUS\n")
        stream.write("========================================\n\n")
        stream.write(f"Complete 17-mode cases used: {len(cases)}.\n")
        stream.write(
            "Axes: Re(lambda)/(2*pi) and Im(lambda)/(2*pi), both in Hz.\n"
        )
        stream.write(
            "The five requested branches are short-period, symmetric bending, "
            "symmetric torsion, antisymmetric bending and antisymmetric torsion.\n"
        )
        stream.write(
            "Elastic signals are reconstructed from the two wing tips after "
            "removing rigid-body heave and pitch; signal derivatives complete "
            "the free-response state used for pole identification.\n"
        )
        stream.write(
            "This is a response-identified closed-loop locus, not a direct "
            "open-loop eigenvalue solution of the MBDyn equations.\n"
        )
        stream.write(
            "A branch in Re(lambda)>0 is unstable. Missing branches were not "
            "sufficiently observable in the imposed excitation.\n"
        )

    print(f"Root locus: {figure_path}")
    print(f"Identified poles: {csv_path}")
    print(f"Method notes: {notes_path}")


if __name__ == "__main__":
    main()
