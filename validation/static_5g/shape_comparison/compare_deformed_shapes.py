#!/usr/bin/env python3
"""Confronta le deformate statiche MBDyn e Nastran per i tre load case."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from cycler import cycler


THESIS_COLORS = [
    "#8B0000",  # dark red
    "#00008B",  # dark blue
    "#66B2FF",  # light blue
    "#006400",  # dark green
    "#CC5500",  # dark orange
]
TEXT_SIZE = 17
TICK_SIZE = 15
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.titlesize": TEXT_SIZE,
    "axes.labelsize": TEXT_SIZE,
    "legend.fontsize": TEXT_SIZE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "axes.prop_cycle": cycler(color=THESIS_COLORS),
})


CASES = {
    1: (
        "case1_symmetric_bending",
        "SYMMETRIC WING BENDING — 10 lbf",
    ),
    2: (
        "case2_antisymmetric_bending",
        "ANTISYMMETRIC WING BENDING — 10 lbf",
    ),
    3: (
        "case3_symmetric_torsion",
        "SYMMETRIC TORSION — 500 lbf-in",
    ),
}
WINGLET_NODES = {
    990021,
    990022,
    990023,
    991021,
    991022,
    991023,
}
FLOAT = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[EDed][+-]?\d+)?"
F06_ROW = re.compile(
    rf"^\s*(\d+)\s+G\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})"
    rf"\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$"
)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Confronta i 45 nodi comuni MBDyn/Nastran e crea grafici PNG."
    )
    parser.add_argument(
        "f06",
        nargs="?",
        type=Path,
        default=root / "nastran/sol101_validation.f06",
        help="F06 prodotto da Nastran",
    )
    parser.add_argument(
        "--mbdyn-dir",
        type=Path,
        default=root / "mbdyn/results",
        help="Cartella contenente i tre file MBDyn .nc",
    )
    parser.add_argument(
        "--nodes",
        type=Path,
        default=root / "nastran/rbe3s.bdf",
        help="File BDF contenente i GRID comuni MBDyn/Nastran",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "plots",
        help="Cartella di destinazione",
    )
    parser.add_argument(
        "--tail-fraction",
        type=float,
        default=0.20,
        help="Frazione finale della storia MBDyn mediata (default: 0.20)",
    )
    return parser.parse_args()


def nastran_float(text: str) -> float:
    """Converte numeri F06 standard; accetta anche l'esponente compatto 1.2-3."""
    value = text.replace("D", "E").replace("d", "E")
    if "E" not in value.upper():
        value = re.sub(r"(?<=\d)([+-]\d+)$", r"E\1", value)
    return float(value)


def load_nodes(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[int, float, float, float]] = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.replace(",", " ").split()
            if len(fields) < 5 or fields[0].upper() != "GRID":
                continue
            node_id = int(fields[1])
            if node_id == 990001 or 990002 <= node_id <= 990023 or 991002 <= node_id <= 991023:
                rows.append(
                    (
                        node_id,
                        nastran_float(fields[-3]),
                        nastran_float(fields[-2]),
                        nastran_float(fields[-1]),
                    )
                )
    if len(rows) != 45:
        raise RuntimeError(
            f"Attesi 45 GRID comuni in {path}, trovati {len(rows)}."
        )
    node_ids = np.array([row[0] for row in rows], dtype=int)
    xyz = np.array([row[1:] for row in rows], dtype=float)
    return node_ids, xyz


def read_f06(path: Path, wanted_nodes: set[int]) -> dict[int, dict[int, np.ndarray]]:
    if not path.is_file():
        raise FileNotFoundError(
            f"F06 non trovato: {path}\n"
            "Eseguire il deck Nastran e copiare qui il relativo .f06."
        )

    results: dict[int, dict[int, np.ndarray]] = {case: {} for case in CASES}
    current_subcase: int | None = None
    in_displacement_table = False
    fatal_messages: list[str] = []

    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "FATAL MESSAGE" in line:
                fatal_messages.append(line.strip())

            subcase_match = re.search(r"\bSUBCASE\s+(\d+)\b", line, re.IGNORECASE)
            if subcase_match:
                current_subcase = int(subcase_match.group(1))

            if "D I S P L A C E M E N T   V E C T O R" in line:
                in_displacement_table = current_subcase in CASES
                continue

            if not in_displacement_table or current_subcase not in CASES:
                continue

            row = F06_ROW.match(line)
            if row:
                node_id = int(row.group(1))
                if node_id in wanted_nodes:
                    results[current_subcase][node_id] = np.array(
                        [nastran_float(row.group(i)) for i in range(2, 8)]
                    )
            elif line.startswith("0") and "PAGE" not in line:
                in_displacement_table = False

    if fatal_messages:
        raise RuntimeError(
            f"Il F06 contiene {len(fatal_messages)} errori fatali Nastran; "
            f"il primo è: {fatal_messages[0]}"
        )

    return results


def read_mbdyn(
    path: Path, node_ids: np.ndarray, tail_fraction: float
) -> tuple[np.ndarray, np.ndarray, float]:
    if not path.is_file():
        raise FileNotFoundError(f"Risultato MBDyn non trovato: {path}")
    if not 0.0 < tail_fraction <= 1.0:
        raise ValueError("--tail-fraction deve essere compreso tra 0 e 1.")

    displacements = np.zeros((len(node_ids), 3))
    rotations = np.zeros((len(node_ids), 3))
    worst_tail_range = 0.0
    with netCDF4.Dataset(path) as dataset:
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        start = max(1, int(np.floor((1.0 - tail_fraction) * len(time))))
        for index, node_id in enumerate(node_ids):
            variable = f"node.struct.{node_id}.X"
            if variable not in dataset.variables:
                raise KeyError(f"Variabile NetCDF assente: {variable}")
            position = np.asarray(dataset.variables[variable][:], dtype=float)
            tail_displacement = position[start:, :3] - position[0, :3]
            displacements[index] = np.mean(tail_displacement, axis=0)
            rotation_variable = f"node.struct.{node_id}.Phi"
            if rotation_variable not in dataset.variables:
                raise KeyError(f"Variabile NetCDF assente: {rotation_variable}")
            rotation = np.asarray(dataset.variables[rotation_variable][:], dtype=float)
            rotations[index] = np.mean(rotation[start:, :3], axis=0)
            worst_tail_range = max(
                worst_tail_range,
                float(np.max(np.ptp(tail_displacement, axis=0))),
            )
    return displacements, rotations, worst_tail_range


def ordered_values(
    node_ids: np.ndarray, xyz: np.ndarray, values: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.lexsort((xyz[:, 0], xyz[:, 2], xyz[:, 1]))
    return node_ids[order], xyz[order], values[order]


def normalized_error_percent(
    mbdyn: np.ndarray, nastran: np.ndarray
) -> tuple[np.ndarray, float, float]:
    """Restituisce errore puntuale normalizzato, RMS ed errore relativo L2."""
    error = mbdyn - nastran
    reference_scale = float(np.max(np.abs(nastran)))
    error_percent = (
        100.0 * error / reference_scale
        if reference_scale > np.finfo(float).eps
        else np.zeros_like(error)
    )
    reference_norm = float(np.linalg.norm(nastran))
    relative_error = (
        float(np.linalg.norm(error) / reference_norm)
        if reference_norm > np.finfo(float).eps
        else float("nan")
    )
    rms = float(np.sqrt(np.mean(error**2)))
    return error_percent, rms, relative_error


def plot_case(
    path: Path,
    title: str,
    node_ids: np.ndarray,
    xyz: np.ndarray,
    mbdyn_displacements: np.ndarray,
    mbdyn_rotations: np.ndarray,
    nastran: np.ndarray,
) -> dict[str, tuple[float, float]]:
    """Crea deformazione in-plane, out-of-plane e torsione con relativi errori."""
    _, ordered_xyz, ordered_mbdyn = ordered_values(
        node_ids, xyz, mbdyn_displacements
    )
    _, _, ordered_mbdyn_rotations = ordered_values(
        node_ids, xyz, mbdyn_rotations
    )
    _, _, ordered_nastran = ordered_values(node_ids, xyz, nastran)
    span = ordered_xyz[:, 1]
    quantities = (
        (
            "in-plane",
            ordered_mbdyn[:, 0],
            ordered_nastran[:, 0],
            "$U_x$ [in]",
        ),
        (
            "out-of-plane",
            ordered_mbdyn[:, 2],
            ordered_nastran[:, 2],
            "$U_z$ [in]",
        ),
        (
            "twist",
            ordered_mbdyn_rotations[:, 1],
            ordered_nastran[:, 4],
            r"$R_y\;[\mathrm{rad}]$",
        ),
    )

    fig, axes = plt.subplots(3, 2, figsize=(15, 13), sharex="col")
    fig.subplots_adjust(
        left=0.08, right=0.98, bottom=0.07, top=0.93, hspace=0.25, wspace=0.20
    )
    metrics: dict[str, tuple[float, float]] = {}
    absolute_data: list[
        tuple[str, np.ndarray, np.ndarray, np.ndarray, str, str]
    ] = []
    for row, (name, mbdyn_values, nastran_values, ylabel) in enumerate(quantities):
        error_percent, rms, relative_error = normalized_error_percent(
            mbdyn_values, nastran_values
        )
        metrics[name] = (rms, relative_error)

        axes[row, 0].plot(
            span, mbdyn_values, color=THESIS_COLORS[1], marker="o",
            markersize=5.5, linewidth=1.8, linestyle="-", label="MBDyn"
        )
        axes[row, 0].plot(
            span, nastran_values, color=THESIS_COLORS[0], marker="s",
            markersize=5.2, linewidth=1.8, linestyle="-", label="Nastran"
        )
        axes[row, 0].axhline(
            0.0, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8
        )
        axes[row, 0].set_ylabel(ylabel)
        axes[row, 0].grid(False)
        axes[row, 0].yaxis.grid(
            True, color="0.82", linewidth=0.7, alpha=0.65
        )
        axes[row, 0].legend(
            frameon=True, fancybox=False, edgecolor="0.25"
        )

        axes[row, 1].plot(
            span, error_percent, color=THESIS_COLORS[2], marker="o",
            markersize=5.5, linewidth=1.8, linestyle="-"
        )
        axes[row, 1].axhline(
            0.0, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8
        )
        axes[row, 1].set_ylabel("Normalized error [%]")
        axes[row, 1].grid(False)
        axes[row, 1].yaxis.grid(
            True, color="0.82", linewidth=0.7, alpha=0.65
        )
        for column in range(2):
            axes[row, column].set_xlabel("$b$ [in]")
            axes[row, column].tick_params(axis="both", labelsize=TICK_SIZE)
            axes[row, column].set_axisbelow(True)

        absolute_error = np.abs(mbdyn_values - nastran_values)
        absolute_unit = "rad" if name == "twist" else "in"
        error_symbol = {
            "in-plane": r"|U_{x,\mathrm{Nast}}-U_{x,\mathrm{MBDyn}}|",
            "out-of-plane": r"|U_{z,\mathrm{Nast}}-U_{z,\mathrm{MBDyn}}|",
            "twist": r"|R_{y,\mathrm{Nast}}-R_{y,\mathrm{MBDyn}}|",
        }[name]
        absolute_data.append(
            (
                ylabel,
                mbdyn_values,
                nastran_values,
                absolute_error,
                absolute_unit,
                error_symbol,
            )
        )

    axes[0, 0].set_title("MBDyn–Nastran comparison")
    axes[0, 1].set_title("Normalized pointwise error")
    fig.suptitle(title, fontsize=TEXT_SIZE, fontweight="normal")
    fig.savefig(path, dpi=450, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    absolute_path = path.with_name(f"{path.stem}_absolute_error{path.suffix}")
    fig, axes = plt.subplots(3, 2, figsize=(15, 13), sharex="col")
    fig.subplots_adjust(
        left=0.08, right=0.98, bottom=0.07, top=0.93, hspace=0.25, wspace=0.20
    )
    for row, (
        ylabel,
        mbdyn_values,
        nastran_values,
        absolute_error,
        absolute_unit,
        error_symbol,
    ) in enumerate(absolute_data):
        axes[row, 0].plot(
            span, mbdyn_values, color=THESIS_COLORS[1], marker="o",
            markersize=5.5, linewidth=1.8, linestyle="-", label="MBDyn"
        )
        axes[row, 0].plot(
            span, nastran_values, color=THESIS_COLORS[0], marker="s",
            markersize=5.2, linewidth=1.8, linestyle="-", label="Nastran"
        )
        axes[row, 0].axhline(
            0.0, color="0.25", linewidth=1.2, linestyle="-", alpha=0.8
        )
        axes[row, 0].set_ylabel(ylabel, fontsize=TEXT_SIZE)
        axes[row, 0].grid(False)
        axes[row, 0].yaxis.grid(
            True, color="0.82", linewidth=0.7, alpha=0.65
        )
        axes[row, 0].legend(
            fontsize=TEXT_SIZE, frameon=True, fancybox=False, edgecolor="0.25"
        )
        axes[row, 0].set_title(
            "MBDyn–Nastran deformation", fontsize=TEXT_SIZE,
            fontweight="normal"
        )

        axes[row, 1].plot(
            span, absolute_error, color=THESIS_COLORS[2], marker="o",
            markersize=5.5, linewidth=1.8, linestyle="-"
        )
        error_max = float(np.max(absolute_error))
        axes[row, 1].set_ylim(0.0, 1.08 * error_max if error_max > 0.0 else 1.0)
        axes[row, 1].margins(x=0.02)
        axes[row, 1].ticklabel_format(
            axis="y", style="sci", scilimits=(-3, 3), useMathText=True
        )
        axes[row, 1].set_ylabel(
            rf"${error_symbol}\;[\mathrm{{{absolute_unit}}}]$",
            fontsize=TEXT_SIZE,
        )
        axes[row, 1].grid(False)
        axes[row, 1].yaxis.grid(
            True, color="0.82", linewidth=0.7, alpha=0.65
        )
        axes[row, 1].set_title(
            "Absolute difference", fontsize=TEXT_SIZE, fontweight="normal"
        )
        for column in range(2):
            axes[row, column].tick_params(axis="both", labelsize=TICK_SIZE)
            axes[row, column].set_xlabel("$b$ [in]", fontsize=TEXT_SIZE)
            axes[row, column].set_axisbelow(True)

    fig.suptitle(title, fontsize=TEXT_SIZE, fontweight="normal")
    fig.savefig(absolute_path, dpi=450, bbox_inches="tight")
    fig.savefig(absolute_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return metrics


def main() -> int:
    args = parse_args()
    node_ids, xyz = load_nodes(args.nodes)
    wanted = set(node_ids.tolist())
    nastran_cases = read_f06(args.f06, wanted)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    keep = np.array([node_id not in WINGLET_NODES for node_id in node_ids])
    comparison_ids = node_ids[keep]
    comparison_xyz = xyz[keep]

    for case_id, (stem, title) in CASES.items():
        missing = wanted.difference(nastran_cases[case_id])
        if missing:
            preview = ", ".join(str(value) for value in sorted(missing)[:8])
            raise RuntimeError(
                f"SUBCASE {case_id}: mancano {len(missing)} nodi nel F06 "
                f"(primi: {preview}). Verificare DISPLACEMENT(PRINT)=101."
            )

        nastran = np.vstack([nastran_cases[case_id][node_id] for node_id in node_ids])
        mbdyn_path = args.mbdyn_dir / f"{stem}.nc"
        mbdyn, mbdyn_rotations, tail_range = read_mbdyn(
            mbdyn_path, node_ids, args.tail_fraction
        )
        metrics = plot_case(
            args.output_dir / f"{stem}_comparison.png",
            title,
            comparison_ids,
            comparison_xyz,
            mbdyn[keep],
            mbdyn_rotations[keep],
            nastran[keep],
        )
        metric_text = ", ".join(
            f"{name}: RMS={rms:.6g} {'deg' if name == 'twist' else 'in'}, "
            f"relative error={relative_error:.3%}"
            for name, (rms, relative_error) in metrics.items()
        )
        print(
            f"SUBCASE {case_id}: {metric_text}, MBDyn tail={tail_range:.3e} in"
        )

    print(f"Confronto completato. Risultati in: {args.output_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileNotFoundError, KeyError, RuntimeError, ValueError) as error:
        print(f"ERRORE: {error}", file=sys.stderr)
        sys.exit(2)
