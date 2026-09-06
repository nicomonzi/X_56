#!/usr/bin/env python3
"""Viewer 3D interattivo delle forme modali dell'X-56.

Comandi:
  mouse sinistro        ruota
  rotella               zoom
  cursore in basso      seleziona direttamente il modo
  pulsanti < e >        modo precedente/successivo
  pulsanti colore       aprono il selettore colori
  frecce destra/sinistra o n/p cambiano modo
  q / Esc               chiude
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

from plot_x56_modes import (
    DEFAULT_BDF,
    DEFAULT_OP2,
    eigenvector_result,
    element_edges,
    mode_displacements,
    node_positions,
    read_x56_bdf,
    result_modes_and_frequencies,
)
from pyNastran.op2.op2 import read_op2


UNDEFORMED_COLOR = "#333333"  # grigio antracite, quasi nero
DEFORMED_COLOR = "#003b7a"    # blu scuro


class ModalViewer:
    def __init__(
        self,
        bdf,
        result,
        modes: np.ndarray,
        frequencies: np.ndarray,
        selected_modes: list[int],
        initial_mode: int,
        positions: dict[int, np.ndarray],
        edges: list[tuple[int, int]],
        deformation: float,
    ) -> None:
        self.bdf = bdf
        self.result = result
        self.modes = modes
        self.frequencies = frequencies
        self.selected_modes = selected_modes
        self.positions = positions
        self.edges = edges
        self.deformation = deformation
        self.mode_to_index = {int(mode): i for i, mode in enumerate(modes)}
        self.current = selected_modes.index(initial_mode)
        self.undeformed_color = UNDEFORMED_COLOR
        self.deformed_color = DEFORMED_COLOR

        self.xyz = np.vstack(list(positions.values()))
        self.span = float(np.max(np.ptp(self.xyz, axis=0)))
        self.undeformed_segments = np.asarray(
            [[positions[n1], positions[n2]] for n1, n2 in edges], dtype=float
        )

        self.fig = plt.figure(figsize=(12, 9), facecolor="white")
        self.fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.17)
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.ax.set_facecolor("white")
        self.ax.set_axis_off()

        self.undeformed_collection = Line3DCollection(
            self.undeformed_segments,
            colors=UNDEFORMED_COLOR,
            linewidths=0.50,
            alpha=0.72,
        )
        self.deformed_collection = Line3DCollection(
            self.undeformed_segments,
            colors=DEFORMED_COLOR,
            linewidths=0.80,
            alpha=1.0,
        )
        self.ax.add_collection3d(self.undeformed_collection)
        self.ax.add_collection3d(self.deformed_collection)

        self._draw_origin_triad()
        self._build_controls()
        self.ax.view_init(elev=24.0, azim=-132.0)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self.fig.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.update_mode()

    def _build_controls(self) -> None:
        slider_ax = self.fig.add_axes((0.23, 0.085, 0.54, 0.035))
        self.mode_slider = Slider(
            slider_ax,
            "Mode",
            valmin=min(self.selected_modes),
            valmax=max(self.selected_modes),
            valinit=self.selected_modes[self.current],
            valstep=np.asarray(self.selected_modes),
            valfmt="%0.0f",
            color=DEFORMED_COLOR,
        )
        self.mode_slider.on_changed(self._on_slider)

        previous_ax = self.fig.add_axes((0.12, 0.075, 0.055, 0.05))
        next_ax = self.fig.add_axes((0.825, 0.075, 0.055, 0.05))
        self.previous_button = Button(previous_ax, "<")
        self.next_button = Button(next_ax, ">")
        self.previous_button.on_clicked(lambda _event: self._step_mode(-1))
        self.next_button.on_clicked(lambda _event: self._step_mode(1))

        undeformed_ax = self.fig.add_axes((0.27, 0.015, 0.21, 0.045))
        deformed_ax = self.fig.add_axes((0.52, 0.015, 0.21, 0.045))
        self.undeformed_color_button = Button(
            undeformed_ax,
            "Colore indeformata",
            color=self.undeformed_color,
            hovercolor="#555555",
        )
        self.deformed_color_button = Button(
            deformed_ax,
            "Colore deformata",
            color=self.deformed_color,
            hovercolor="#15559a",
        )
        for button in (
            self.undeformed_color_button,
            self.deformed_color_button,
        ):
            button.label.set_color("white")
        self.undeformed_color_button.on_clicked(
            lambda _event: self._choose_color("undeformed")
        )
        self.deformed_color_button.on_clicked(
            lambda _event: self._choose_color("deformed")
        )

    def _choose_color(self, target: str) -> None:
        """Apre il color picker nativo e applica subito il colore scelto."""
        try:
            import tkinter as tk
            from tkinter import colorchooser
        except ImportError:
            print("Color picker non disponibile: manca tkinter.")
            return

        current = (
            self.undeformed_color if target == "undeformed"
            else self.deformed_color
        )
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        _rgb, selected = colorchooser.askcolor(
            color=current, parent=root, title="Scegli il colore"
        )
        root.destroy()
        if selected is None:
            return

        if target == "undeformed":
            self.undeformed_color = selected
            self.undeformed_collection.set_color(selected)
            self.undeformed_color_button.ax.set_facecolor(selected)
        else:
            self.deformed_color = selected
            self.deformed_collection.set_color(selected)
            self.deformed_color_button.ax.set_facecolor(selected)
            self.mode_slider.poly.set_color(selected)
        self.fig.canvas.draw_idle()

    def _draw_origin_triad(self) -> None:
        """Disegna la terna nel punto (0,0,0) del sistema basico Nastran."""
        length = 0.12 * self.span
        origin = np.zeros(3)
        for axis, (color, label) in enumerate(
            zip(("#d62728", "#2ca02c", "#1565c0"), ("X", "Y", "Z"))
        ):
            direction = np.zeros(3)
            direction[axis] = length
            self.ax.quiver(
                *origin,
                *direction,
                color=color,
                linewidth=2.0,
                arrow_length_ratio=0.14,
            )
            endpoint = origin + 1.10 * direction
            self.ax.text(
                *endpoint,
                label,
                color=color,
                fontsize=11,
                ha="center",
                va="center",
            )

    def _set_equal_limits(self, xyz: np.ndarray) -> None:
        xyz_min = xyz.min(axis=0)
        xyz_max = xyz.max(axis=0)
        center = 0.5 * (xyz_min + xyz_max)
        half_range = 0.54 * float(np.max(xyz_max - xyz_min))
        self.ax.set_xlim(center[0] - half_range, center[0] + half_range)
        self.ax.set_ylim(center[1] - half_range, center[1] + half_range)
        self.ax.set_zlim(center[2] - half_range, center[2] + half_range)
        self.ax.set_box_aspect((1.0, 1.0, 1.0))

    def update_mode(self) -> None:
        mode = self.selected_modes[self.current]
        result_index = self.mode_to_index[mode]
        displacements = mode_displacements(self.result, self.bdf, result_index)
        common = set(self.positions).intersection(displacements)
        max_displacement = max(
            np.linalg.norm(displacements[nid]) for nid in common
        )
        scale = self.deformation * self.span / max_displacement

        deformed = {
            nid: position + scale * displacements.get(nid, np.zeros(3))
            for nid, position in self.positions.items()
        }
        segments = np.asarray(
            [[deformed[n1], deformed[n2]] for n1, n2 in self.edges], dtype=float
        )
        self.deformed_collection.set_segments(segments)

        all_xyz = np.vstack(
            (self.undeformed_segments.reshape(-1, 3), segments.reshape(-1, 3))
        )
        self._set_equal_limits(all_xyz)
        frequency = float(self.frequencies[result_index])
        self.ax.set_title(
            f"Mode {mode} — {frequency:.3f} Hz", fontsize=15, pad=12
        )
        self.fig.canvas.draw_idle()
        print(f"Mode {mode} — {frequency:.3f} Hz")

    def _on_slider(self, value: float) -> None:
        mode = int(round(value))
        self.current = self.selected_modes.index(mode)
        self.update_mode()

    def _step_mode(self, step: int) -> None:
        self.current = (self.current + step) % len(self.selected_modes)
        self.mode_slider.set_val(self.selected_modes[self.current])

    def _on_key(self, event) -> None:
        if event.key in ("right", "n"):
            self._step_mode(1)
        elif event.key in ("left", "p"):
            self._step_mode(-1)
        elif event.key in ("q", "escape"):
            plt.close(self.fig)

    def _on_scroll(self, event) -> None:
        """Zoom centrato sulla vista, indipendente dalla toolbar grafica."""
        if event.inaxes is not self.ax:
            return
        factor = 0.88 if event.button == "up" else 1.14
        for getter, setter in (
            (self.ax.get_xlim3d, self.ax.set_xlim3d),
            (self.ax.get_ylim3d, self.ax.set_ylim3d),
            (self.ax.get_zlim3d, self.ax.set_zlim3d),
        ):
            low, high = getter()
            center = 0.5 * (low + high)
            half = 0.5 * (high - low) * factor
            setter(center - half, center + half)
        self.fig.canvas.draw_idle()

    def show(self) -> None:
        plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Viewer 3D interattivo delle forme modali X-56."
    )
    parser.add_argument("--bdf", type=Path, default=DEFAULT_BDF)
    parser.add_argument("--op2", type=Path, default=DEFAULT_OP2)
    parser.add_argument(
        "--first-mode",
        type=int,
        default=1,
        help="Modo mostrato all'apertura (default: 1, modo rigido).",
    )
    parser.add_argument(
        "--count",
        type=int,
        default=31,
        help="Numero di modi caricati a partire dal primo (default: tutti i 31).",
    )
    parser.add_argument(
        "--deformation",
        type=float,
        default=0.12,
        help="Ampiezza grafica della deformazione (default: 0.12).",
    )
    parser.add_argument(
        "--wireframe-stride",
        type=int,
        default=3,
        help="Alleggerimento linee interne della skin: 1=mesh completa, "
             "3=default più leggero.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.bdf.is_file():
        raise SystemExit(f"BDF non trovato: {args.bdf}")
    if not args.op2.is_file():
        raise SystemExit(f"OP2 non trovato: {args.op2}")
    if args.count < 1 or args.deformation <= 0.0 or args.wireframe_stride < 1:
        raise SystemExit(
            "--count, --deformation e --wireframe-stride devono essere positivi."
        )

    print(f"Leggo geometria: {args.bdf}")
    bdf = read_x56_bdf(args.bdf)
    print(f"Leggo risultati: {args.op2}")
    op2 = read_op2(
        str(args.op2), combine=True, build_dataframe=False, debug=False
    )
    result = eigenvector_result(op2)
    modes, frequencies = result_modes_and_frequencies(result)

    if args.first_mode not in modes:
        raise SystemExit(
            f"Il modo {args.first_mode} non è presente; disponibili: {modes.tolist()}"
        )
    selected = modes[:args.count].astype(int).tolist()
    if args.first_mode not in selected:
        selected.append(args.first_mode)
        selected.sort()

    positions = node_positions(bdf)
    # Solo elementi shell: niente RBE, barre, molle o collegamenti ausiliari.
    # Il risultato è una skin a wireframe molto più leggera da ruotare.
    edges = element_edges(
        bdf,
        set(positions),
        shell_only=True,
        wireframe_stride=args.wireframe_stride,
    )
    print(f"Skin wireframe: {len(edges)} spigoli shell")
    viewer = ModalViewer(
        bdf=bdf,
        result=result,
        modes=modes,
        frequencies=frequencies,
        selected_modes=selected,
        initial_mode=args.first_mode,
        positions=positions,
        edges=edges,
        deformation=args.deformation,
    )
    viewer.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
