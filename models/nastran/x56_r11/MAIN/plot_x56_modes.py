#!/usr/bin/env python3
"""Plot delle forme modali X-56: geometria indeformata e deformata in 3D.

Legge la geometria dal BDF (inclusi tutti gli INCLUDE) e gli autovettori
direttamente dall'OP2. Per default disegna i primi 10 modi elastici.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

try:
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.op2 import read_op2
except ImportError:
    sys.exit(
        "Manca pyNastran. Installalo con:\n"
        "  python3 -m pip install pyNastran\n"
        "Poi rilancia questo script."
    )


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BDF = SCRIPT_DIR / "nsvibe_cfg24611_10lb.bdf"
DEFAULT_OP2 = SCRIPT_DIR / "nsvibe_cfg24611_10lb.op2"

# Elementi utili per mostrare la forma della struttura. Gli elementi rigidi
# sono inclusi perché rendono leggibili collegamenti e superfici di comando.
LINE_TYPES = {
    "CROD", "CONROD", "CTUBE", "CBAR", "CBEAM", "CBUSH", "CBUSH1D",
    "CELAS1", "CELAS2", "CELAS3", "CELAS4", "RBE1", "RBE2", "RBE3",
}
SHELL_TYPES = {
    "CTRIA3", "CTRIA6", "CQUAD4", "CQUAD8", "CQUADR", "CTRIAR", "CSHEAR",
}


def parse_mode_spec(spec: str) -> list[int]:
    """Converte, per esempio, '7-12,15,18' in una lista di modi."""
    modes: list[int] = []
    for item in spec.replace(" ", "").split(","):
        if not item:
            continue
        if "-" in item:
            first, last = (int(value) for value in item.split("-", 1))
            step = 1 if last >= first else -1
            modes.extend(range(first, last + step, step))
        else:
            modes.append(int(item))
    return list(dict.fromkeys(modes))


def eigenvector_result(op2):
    """Restituisce l'unico set di autovettori, con errore chiaro se ambiguo."""
    if not op2.eigenvectors:
        raise RuntimeError("L'OP2 non contiene autovettori (tabella OUG).")
    if len(op2.eigenvectors) > 1:
        keys = list(op2.eigenvectors)
        print(f"Più subcase trovati {keys}: uso il primo ({keys[0]}).")
    return next(iter(op2.eigenvectors.values()))


def read_x56_bdf(filename: Path) -> BDF:
    """Legge il deck NASA, tollerando gli ID proprietà duplicati presenti."""
    bdf = BDF(debug=False)
    # Il deck originale contiene almeno un ID condiviso da PBAR e PSHELL.
    # MSC Nastran lo accetta; pyNastran lo segnala come errore recuperabile.
    bdf.set_error_storage(
        nparse_errors=1000,
        stop_on_parsing_error=False,
        nxref_errors=1000,
        stop_on_xref_error=False,
    )
    bdf.read_bdf(str(filename), validate=False, xref=True)
    return bdf


def result_modes_and_frequencies(result) -> tuple[np.ndarray, np.ndarray]:
    modes = np.asarray(result.modes, dtype=int)
    for name in ("mode_cycles", "cycles", "freqs"):
        values = getattr(result, name, None)
        if values is not None and len(values) == len(modes):
            return modes, np.asarray(values, dtype=float)
    eigenvalues = np.asarray(result.eigns, dtype=float)
    frequencies = np.sqrt(np.abs(eigenvalues)) / (2.0 * np.pi)
    return modes, frequencies


def automatic_modes(modes: np.ndarray, frequencies: np.ndarray, count: int) -> list[int]:
    # I modi quasi-rigidi del modello free-free sono prossimi a 0 Hz.
    elastic = modes[frequencies >= 0.1]
    return elastic[:count].astype(int).tolist()


def node_positions(bdf) -> dict[int, np.ndarray]:
    positions: dict[int, np.ndarray] = {}
    for nid, node in bdf.nodes.items():
        try:
            positions[nid] = np.asarray(node.get_position(), dtype=float)
        except Exception:
            # SPOINT/EPOINT non possiedono una posizione 3D.
            continue
    return positions


def element_edges(
    bdf,
    available_nodes: set[int],
    shell_only: bool = False,
    wireframe_stride: int = 1,
) -> list[tuple[int, int]]:
    """Estrae gli spigoli della mesh senza duplicati."""
    edge_counts: Counter[tuple[int, int]] = Counter()
    for element in bdf.elements.values():
        etype = element.type
        if etype not in SHELL_TYPES and (shell_only or etype not in LINE_TYPES):
            continue
        try:
            nids = [nid for nid in element.node_ids if nid in available_nodes]
        except Exception:
            continue
        if len(nids) < 2:
            continue

        if etype in SHELL_TYPES:
            # Per elementi quadratici conserva anche i nodi intermedi.
            pairs = zip(nids, nids[1:] + nids[:1])
        elif etype in {"RBE1", "RBE2", "RBE3"}:
            # Una stella dal nodo di riferimento evita una polilinea fittizia.
            pairs = ((nids[0], nid) for nid in nids[1:])
        else:
            pairs = ((nids[0], nids[1]),)

        for n1, n2 in pairs:
            if n1 != n2:
                edge_counts[(min(n1, n2), max(n1, n2))] += 1
    if not edge_counts:
        kind = "shell" if shell_only else "plottabile"
        raise RuntimeError(f"Non è stato trovato alcun elemento {kind} nel BDF.")

    edges = sorted(edge_counts)
    if shell_only and wireframe_stride > 1:
        # Mantiene tutti i bordi liberi della skin e alleggerisce solo le linee
        # interne della mesh. La sagoma rimane completa, ma il rendering 3D è
        # sensibilmente più rapido.
        boundary = [edge for edge in edges if edge_counts[edge] == 1]
        interior = [edge for edge in edges if edge_counts[edge] > 1]
        edges = sorted(set(boundary).union(interior[::wireframe_stride]))
    return edges


def mode_displacements(result, bdf, time_index: int) -> dict[int, np.ndarray]:
    """Traslazioni dell'autovettore espresse nel sistema basico Nastran."""
    data = np.real(np.asarray(result.data[time_index, :, :3], dtype=float))
    node_ids = np.asarray(result.node_gridtype[:, 0], dtype=int)
    displacements: dict[int, np.ndarray] = {}

    for nid, vector in zip(node_ids, data):
        node = bdf.nodes.get(int(nid))
        if node is None or node.type != "GRID":
            continue
        # Gli output T1,T2,T3 sono espressi nel CD del GRID. Occorre ruotarli
        # nel sistema basico prima di sommarli alle coordinate globali.
        try:
            cd_ref = node.cd_ref
            if node.Cd() != 0 and cd_ref is not None:
                vector = cd_ref.transform_vector_to_global(vector)
        except (AttributeError, TypeError):
            pass
        displacements[int(nid)] = np.asarray(vector, dtype=float)
    return displacements


def equal_3d_axes(ax, xyz: np.ndarray, padding: float = 0.06) -> None:
    xyz_min = xyz.min(axis=0)
    xyz_max = xyz.max(axis=0)
    center = 0.5 * (xyz_min + xyz_max)
    half_range = 0.5 * float(np.max(xyz_max - xyz_min)) * (1.0 + 2.0 * padding)
    if half_range == 0.0:
        half_range = 1.0
    ax.set_xlim(center[0] - half_range, center[0] + half_range)
    ax.set_ylim(center[1] - half_range, center[1] + half_range)
    ax.set_zlim(center[2] - half_range, center[2] + half_range)
    ax.set_box_aspect((1.0, 1.0, 1.0))


def add_reference_triad(ax, xyz: np.ndarray) -> None:
    xyz_min = xyz.min(axis=0)
    xyz_max = xyz.max(axis=0)
    length = 0.12 * float(np.max(xyz_max - xyz_min))
    origin = xyz_min + 0.04 * (xyz_max - xyz_min)
    colors = ("#d62728", "#2ca02c", "#1f77b4")
    labels = ("X", "Y", "Z")
    for axis, (color, label) in enumerate(zip(colors, labels)):
        direction = np.zeros(3)
        direction[axis] = length
        ax.quiver(*origin, *direction, color=color, linewidth=1.8,
                  arrow_length_ratio=0.14)
        end = origin + 1.10 * direction
        ax.text(*end, label, color=color, fontsize=10, ha="center", va="center")


def clean_white_axes(ax) -> None:
    ax.set_facecolor("white")
    ax.figure.set_facecolor("white")
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 1.0))
        axis.pane.set_edgecolor((0.78, 0.78, 0.78, 1.0))
        axis._axinfo["grid"]["color"] = (0.86, 0.86, 0.86, 1.0)
        axis._axinfo["grid"]["linewidth"] = 0.55
    ax.grid(True)
    ax.set_xlabel("X [in]")
    ax.set_ylabel("Y [in]")
    ax.set_zlabel("Z [in]")


def plot_mode(
    mode: int,
    frequency: float,
    positions: dict[int, np.ndarray],
    edges: list[tuple[int, int]],
    displacements: dict[int, np.ndarray],
    scale_fraction: float,
    output_dir: Path,
    dpi: int,
) -> Path:
    common = set(positions).intersection(displacements)
    if not common:
        raise RuntimeError(f"Nessun GRID comune tra geometria e modo {mode}.")

    span = float(np.max(np.ptp(np.vstack(list(positions.values())), axis=0)))
    max_modal_displacement = max(np.linalg.norm(displacements[nid]) for nid in common)
    scale = scale_fraction * span / max_modal_displacement

    deformed = {
        nid: xyz + scale * displacements.get(nid, np.zeros(3))
        for nid, xyz in positions.items()
    }
    undeformed_segments = np.asarray(
        [[positions[n1], positions[n2]] for n1, n2 in edges], dtype=float
    )
    deformed_segments = np.asarray(
        [[deformed[n1], deformed[n2]] for n1, n2 in edges], dtype=float
    )

    fig = plt.figure(figsize=(12, 9), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")
    clean_white_axes(ax)
    ax.add_collection3d(
        Line3DCollection(
            undeformed_segments, colors="#9a9a9a", linewidths=0.45, alpha=0.80
        )
    )
    ax.add_collection3d(
        Line3DCollection(
            deformed_segments, colors="#0755c9", linewidths=0.70, alpha=0.98
        )
    )

    all_xyz = np.vstack((undeformed_segments.reshape(-1, 3),
                         deformed_segments.reshape(-1, 3)))
    equal_3d_axes(ax, all_xyz)
    add_reference_triad(ax, all_xyz)
    ax.view_init(elev=24.0, azim=-132.0)
    ax.set_title(f"Mode {mode} — {frequency:.3f} Hz", fontsize=15, pad=18)

    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"mode_{mode:02d}_{frequency:.3f}Hz.png"
    fig.savefig(output, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot 3D delle forme modali dell'X-56 da BDF + OP2."
    )
    parser.add_argument("--bdf", type=Path, default=DEFAULT_BDF)
    parser.add_argument("--op2", type=Path, default=DEFAULT_OP2)
    parser.add_argument(
        "--modes",
        help="Modi da plottare, es. '7-16' o '7,9,12'. "
             "Default: primi modi elastici.",
    )
    parser.add_argument(
        "--count", type=int, default=10,
        help="Numero di modi elastici automatici (default: 10).",
    )
    parser.add_argument(
        "--deformation", type=float, default=0.12,
        help="Spostamento massimo come frazione della dimensione del modello "
             "(default: 0.12).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=SCRIPT_DIR / "modal_plots_10lb"
    )
    parser.add_argument("--dpi", type=int, default=220)
    args = parser.parse_args()

    if not args.bdf.is_file():
        parser.error(f"BDF non trovato: {args.bdf}")
    if not args.op2.is_file():
        parser.error(f"OP2 non trovato: {args.op2}")
    if args.count < 1 or args.deformation <= 0.0 or args.dpi < 50:
        parser.error("--count deve essere positivo, --deformation > 0 e --dpi >= 50.")

    print(f"Leggo geometria: {args.bdf}")
    bdf = read_x56_bdf(args.bdf)
    print(f"Leggo risultati: {args.op2}")
    op2 = read_op2(str(args.op2), combine=True, build_dataframe=False, debug=False)
    result = eigenvector_result(op2)
    modes, frequencies = result_modes_and_frequencies(result)

    requested = (
        parse_mode_spec(args.modes)
        if args.modes
        else automatic_modes(modes, frequencies, args.count)
    )
    available = set(modes.tolist())
    missing = [mode for mode in requested if mode not in available]
    if missing:
        parser.error(f"Modi non presenti nell'OP2: {missing}; disponibili: {modes.tolist()}")

    positions = node_positions(bdf)
    edges = element_edges(bdf, set(positions))
    mode_to_index = {int(mode): index for index, mode in enumerate(modes)}
    print(
        f"{len(positions)} GRID, {len(edges)} spigoli; "
        f"plot dei modi {requested}"
    )

    for mode in requested:
        index = mode_to_index[mode]
        displacements = mode_displacements(result, bdf, index)
        output = plot_mode(
            mode=mode,
            frequency=float(frequencies[index]),
            positions=positions,
            edges=edges,
            displacements=displacements,
            scale_fraction=args.deformation,
            output_dir=args.output_dir,
            dpi=args.dpi,
        )
        print(f"Creato: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
