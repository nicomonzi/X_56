#!/usr/bin/env python3
"""Interactive 3D modal-shape viewer; images are saved only on user request."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from compare_mode_shapes import SPAN_NODE_IDS, parse_fem


class InteractiveModalPlot:
    """Display and interactively arrange one modal shape."""

    def __init__(
        self,
        mode: int,
        frequency: float,
        coordinates: np.ndarray,
        shape: np.ndarray,
        amplitude: float,
        centre_at_root: bool,
        screenshot_dir: Path,
    ) -> None:
        self.mode = mode
        self.frequency = frequency
        self.coordinates = coordinates
        self.amplitude = amplitude
        self.screenshot_dir = screenshot_dir
        self.zoom = 0.92
        self.linewidth = 2
        self.markersize = 2.4
        self.axes_position = np.array([0.08, 0.10, 0.86, 0.80])

        display_shape = shape[:, :3].copy()
        if centre_at_root:
            # Change only the displayed reference frame: the relative modal
            # deformation and the source modal data remain unchanged.
            display_shape -= display_shape[0]
        scale = float(np.max(np.linalg.norm(display_shape, axis=1)))
        if scale <= 0.0:
            raise ValueError(f"Mode {mode} has no translational displacement")
        display_shape /= scale
        self.deformed = coordinates + amplitude * display_shape

        self.undeformed_right = np.vstack((coordinates[0], coordinates[23:]))
        self.deformed_right = np.vstack((self.deformed[0], self.deformed[23:]))
        displayed_points = np.vstack((coordinates, self.deformed))
        lower = displayed_points.min(axis=0) - 10.0
        upper = displayed_points.max(axis=0) + 10.0
        self.plot_range = upper - lower
        self.limits = tuple(zip(lower, upper))

        self.figure = plt.figure(figsize=(9.0, 5.2))
        self.axis = self.figure.add_subplot(111, projection="3d")
        self.axis.set_position(self.axes_position)
        self.lines = self._draw_lines()
        self._format_axes()
        self._print_help()

        self.figure.canvas.mpl_connect("key_press_event", self._on_key)
        self.figure.canvas.mpl_connect("scroll_event", self._on_scroll)

    def _draw_lines(self):
        line_u_left, = self.axis.plot(
            self.coordinates[:23, 0],
            self.coordinates[:23, 1],
            self.coordinates[:23, 2],
            "o-",
            color="0.40",
            linewidth=self.linewidth,
            markersize=self.markersize,
            label="Undeformed",
        )
        line_u_right, = self.axis.plot(
            self.undeformed_right[:, 0],
            self.undeformed_right[:, 1],
            self.undeformed_right[:, 2],
            "o-",
            color="0.40",
            linewidth=self.linewidth,
            markersize=self.markersize,
        )
        line_d_left, = self.axis.plot(
            self.deformed[:23, 0],
            self.deformed[:23, 1],
            self.deformed[:23, 2],
            "o-",
            color="tab:blue",
            linewidth=self.linewidth,
            markersize=self.markersize,
            label="Deformed",
        )
        line_d_right, = self.axis.plot(
            self.deformed_right[:, 0],
            self.deformed_right[:, 1],
            self.deformed_right[:, 2],
            "o-",
            color="tab:blue",
            linewidth=self.linewidth,
            markersize=self.markersize,
        )
        return line_u_left, line_u_right, line_d_left, line_d_right

    def _format_axes(self) -> None:
        self.axis.set_xlabel("X [in]", labelpad=5)
        self.axis.set_ylabel("Y [in]", labelpad=18)
        self.axis.set_zlabel("Z [in]", labelpad=6)
        self.axis.xaxis.set_major_locator(MaxNLocator(nbins=4))
        self.axis.yaxis.set_major_locator(MaxNLocator(nbins=7))
        self.axis.zaxis.set_major_locator(MaxNLocator(nbins=5))
        self.axis.tick_params(axis="both", which="major", labelsize=9, pad=1)
        self.axis.set_title(
            f"Mode {self.mode} — {self.frequency:.4f} Hz",
            fontsize=13,
            pad=4,
        )
        self.axis.set_xlim(*self.limits[0])
        self.axis.set_ylim(*self.limits[1])
        self.axis.set_zlim(*self.limits[2])
        self.axis.set_box_aspect(tuple(self.plot_range), zoom=self.zoom)
        self.axis.set_proj_type("ortho")
        self.axis.view_init(elev=25.0, azim=150.0)

        self.legend = self.axis.legend(
            loc="lower center",
            bbox_to_anchor=(0.58, 0.01),
            ncol=2,
            fontsize=9,
            frameon=False,
        )
        self.legend.set_draggable(True, use_blit=True)

    @staticmethod
    def _print_help() -> None:
        print(
            "\nInteractive controls:\n"
            "  mouse drag        rotate the 3D view\n"
            "  drag legend       move the legend freely\n"
            "  arrow keys        move the complete graph in the canvas\n"
            "  + / - or wheel    zoom in / out\n"
            "  [ / ]             decrease / increase line thickness\n"
            "  , / .             decrease / increase marker size\n"
            "  r                  reset view and graph position\n"
            "  s                  save a PNG screenshot now\n"
            "  Esc                close the viewer\n"
            "The toolbar Save button can also be used to choose a filename.\n"
        )

    def _refresh_style(self) -> None:
        for line in self.lines:
            line.set_linewidth(self.linewidth)
            line.set_markersize(self.markersize)
        self.axis.set_box_aspect(tuple(self.plot_range), zoom=self.zoom)
        self.figure.canvas.draw_idle()

    def _move_axes(self, dx: float, dy: float) -> None:
        self.axes_position[:2] += (dx, dy)
        self.axis.set_position(self.axes_position)
        self.figure.canvas.draw_idle()

    def _save_screenshot(self) -> None:
        self.screenshot_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = self.screenshot_dir / f"mode_{self.mode:02d}_{stamp}.png"
        self.figure.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.05)
        print(f"Screenshot saved: {path}")

    def _reset(self) -> None:
        self.zoom = 0.92
        self.linewidth = 0.85
        self.markersize = 2.4
        self.axes_position = np.array([0.08, 0.10, 0.86, 0.80])
        self.axis.set_position(self.axes_position)
        self.axis.view_init(elev=25.0, azim=150.0)
        self._refresh_style()

    def _on_scroll(self, event) -> None:
        self.zoom = np.clip(
            self.zoom * (1.06 if event.button == "up" else 1.0 / 1.06),
            0.55,
            1.55,
        )
        self._refresh_style()

    def _on_key(self, event) -> None:
        key = event.key
        if key == "left":
            self._move_axes(-0.015, 0.0)
        elif key == "right":
            self._move_axes(0.015, 0.0)
        elif key == "up":
            self._move_axes(0.0, 0.015)
        elif key == "down":
            self._move_axes(0.0, -0.015)
        elif key in {"+", "="}:
            self.zoom = min(1.55, self.zoom * 1.06)
            self._refresh_style()
        elif key in {"-", "_"}:
            self.zoom = max(0.55, self.zoom / 1.06)
            self._refresh_style()
        elif key == "[":
            self.linewidth = max(0.25, self.linewidth - 0.10)
            self._refresh_style()
        elif key == "]":
            self.linewidth = min(3.0, self.linewidth + 0.10)
            self._refresh_style()
        elif key == ",":
            self.markersize = max(0.0, self.markersize - 0.25)
            self._refresh_style()
        elif key == ".":
            self.markersize = min(8.0, self.markersize + 0.25)
            self._refresh_style()
        elif key == "r":
            self._reset()
        elif key == "s":
            self._save_screenshot()
        elif key == "escape":
            plt.close(self.figure)

    def show(self) -> None:
        plt.show()


def main() -> int:
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", type=int, default=7, help="Nastran/FEM mode number")
    parser.add_argument(
        "--fem",
        type=Path,
        default=base.parent / "NASTRAN/FEMGEN40/mbdyn_modal.fem",
    )
    parser.add_argument(
        "--display-amplitude",
        type=float,
        default=25.0,
        help="Maximum displayed displacement in inches",
    )
    parser.add_argument(
        "--no-center",
        action="store_true",
        help="Do not make the central undeformed/deformed nodes coincide",
    )
    parser.add_argument(
        "--screenshot-dir",
        type=Path,
        default=base / "manual_screenshots",
        help="Destination used only when the S key is pressed",
    )
    args = parser.parse_args()

    coordinates, shapes, frequencies = parse_fem(args.fem.resolve(), {args.mode})
    if args.mode not in shapes or args.mode not in frequencies:
        raise ValueError(f"Mode {args.mode} is not available in {args.fem}")
    coordinate_array = np.vstack(
        [coordinates[node_id] for node_id in SPAN_NODE_IDS]
    )
    viewer = InteractiveModalPlot(
        mode=args.mode,
        frequency=frequencies[args.mode],
        coordinates=coordinate_array,
        shape=shapes[args.mode],
        amplitude=args.display_amplitude,
        centre_at_root=not args.no_center,
        screenshot_dir=args.screenshot_dir.resolve(),
    )
    viewer.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
