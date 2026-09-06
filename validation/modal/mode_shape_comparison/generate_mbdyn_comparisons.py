#!/usr/bin/env python3
"""Run MBDyn modes 7-18 and compare MOV reconstructions with the Nastran FEM export."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from compare_mode_shapes import (
    SPAN_NODE_IDS,
    modal_assurance_criterion,
    normalized_and_aligned,
    parse_fem,
    parse_mbdyn_mov,
    parse_modes,
    spanwise_arc_length,
)

MODE_NAMES = {
    7: ("SWB1", "First Symmetric Wing Bending"),
    8: ("AWB1", "First Antisymmetric Wing Bending"),
    9: ("SWT1", "First Symmetric Wing Torsion"),
    10: ("SFA", "Symmetric Fore-Aft"),
    11: ("AWT1", "First Antisymmetric Wing Torsion"),
    12: ("SWB2", "Second Symmetric Wing Bending"),
    13: ("AMLGL", "Antisymmetric Main Landing Gear Lateral"),
    14: ("SWL", "Symmetric Winglet Lateral"),
    15: ("BOOM-H", "Horizontal Boom Bending"),
    16: ("LOCAL-1", "Local Structural Mode"),
    17: ("LOCAL-2", "Local Structural Mode"),
    18: ("AWB2", "Second Antisymmetric Wing Bending"),
}


def main() -> int:
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", default="7-18")
    parser.add_argument("--mbdyn", default="/usr/local/mbdyn/bin/mbdyn")
    parser.add_argument("--output", type=Path, default=base / "comparison_results")
    parser.add_argument(
        "--display-amplitude",
        type=float,
        default=25.0,
        help="Maximum displayed 3D modal displacement in inches",
    )
    args = parser.parse_args()

    modes = parse_modes(args.modes)
    fem_path = (base / "../NASTRAN/FEMGEN40/mbdyn_modal.fem").resolve()
    template_path = base / "mbdyn/modal_mode_check.mbd"
    nodes_path = (base / "mbdyn/span_nodes.mbd").resolve()
    coordinates, fem_shapes, frequencies = parse_fem(fem_path, set(modes))
    span = spanwise_arc_length(coordinates)
    coordinate_array = np.vstack([coordinates[node_id] for node_id in SPAN_NODE_IDS])
    order = np.argsort(span)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    template = template_path.read_text(encoding="utf-8")
    summary = []

    with tempfile.TemporaryDirectory(prefix="x56_modal_comparison_") as temporary:
        temporary_path = Path(temporary)
        for mode in modes:
            case = re.sub(
                r"set: const integer MODE_TO_EXCITE = \d+;",
                f"set: const integer MODE_TO_EXCITE = {mode};",
                template,
            )
            case = case.replace('include: "span_nodes.mbd";', f'include: "{nodes_path}";')
            case = case.replace(
                '"../../NASTRAN/FEMGEN40/mbdyn_modal.fem"', f'"{fem_path}"'
            )
            input_path = temporary_path / f"mode_{mode:02d}.mbd"
            input_path.write_text(case, encoding="utf-8")
            prefix = temporary_path / f"mode_{mode:02d}"
            result = subprocess.run(
                [args.mbdyn, "-f", str(input_path), "-o", str(prefix)],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode:
                raise RuntimeError(
                    f"MBDyn failed for mode {mode}:\n{result.stdout}\n{result.stderr}"
                )

            reconstructed = parse_mbdyn_mov(prefix.with_suffix(".mov"))
            reference, mbdyn = normalized_and_aligned(fem_shapes[mode], reconstructed)
            difference = reference - mbdyn
            mac = modal_assurance_criterion(reference, mbdyn)
            rms = float(np.sqrt(np.mean(difference**2)))
            maximum = float(np.max(np.linalg.norm(difference, axis=1)))
            summary.append([mode, frequencies[mode], mac, rms, maximum])

            fig, axes = plt.subplots(3, 1, figsize=(9.0, 10.0), sharex=True)
            labels = (
                "Normalized X displacement",
                "Normalized Y displacement",
                "Normalized Z displacement",
            )
            for component, (axis, label) in enumerate(zip(axes, labels)):
                axis.plot(
                    span[order],
                    reference[order, component],
                    "o-",
                    label="Nastran modal export (.fem)",
                )
                axis.plot(
                    span[order],
                    mbdyn[order, component],
                    "s--",
                    label="MBDyn modal-joint reconstruction (.mov)",
                )
                axis.plot(
                    span[order],
                    difference[order, component],
                    ":",
                    color="0.35",
                    label="Difference",
                )
                axis.set_ylabel(label)
                axis.grid(True, alpha=0.3)
            axes[0].legend(loc="best", fontsize=8)
            axes[-1].set_xlabel("Signed spanwise arc length [in]")
            fig.suptitle(f"Mode {mode} — {frequencies[mode]:.4f} Hz")
            fig.tight_layout()
            fig.savefig(output / f"mode_{mode:02d}_nastran_vs_mbdyn.png", dpi=180)
            plt.close(fig)

            # MBDyn physical reconstruction in the original X/Y/Z frame.
            # The modal vector is displayed with a common maximum amplitude so
            # modes remain visually comparable.
            # Display the flexible deformation in a translating frame whose
            # origin follows the common centre/root node. This makes the
            # undeformed and deformed centre nodes coincide without changing
            # the numerical modal comparison performed above.
            display_shape = mbdyn - mbdyn[0]
            display_norm = float(
                np.max(np.linalg.norm(display_shape, axis=1))
            )
            if display_norm > 0.0:
                display_shape /= display_norm
            deformed = coordinate_array + args.display_amplitude * display_shape
            undeformed_right = np.vstack((coordinate_array[0], coordinate_array[23:]))
            right_with_root = np.vstack((deformed[0], deformed[23:]))
            fig = plt.figure(figsize=(8.0, 3.5))
            axis_3d = fig.add_subplot(111, projection="3d")
            axis_3d.plot(
                coordinate_array[:23, 0],
                coordinate_array[:23, 1],
                coordinate_array[:23, 2],
                "o-",
                color="0.40",
                linewidth=1.5,
                markersize=3.0,
                label="Undeformed",
            )
            axis_3d.plot(
                undeformed_right[:, 0],
                undeformed_right[:, 1],
                undeformed_right[:, 2],
                "o-",
                color="0.40",
                linewidth=1.5,
                markersize=3.0,
            )
            axis_3d.plot(
                deformed[:23, 0],
                deformed[:23, 1],
                deformed[:23, 2],
                "o-",
                color="tab:blue",
                linewidth=2.0,
                markersize=4.0,
                label="Deformed",
            )
            axis_3d.plot(
                right_with_root[:, 0],
                right_with_root[:, 1],
                right_with_root[:, 2],
                "o-",
                color="tab:blue",
                linewidth=2.0,
                markersize=4.0,
            )
            # A fixed 2D placement keeps the short, receding X-axis label
            # inside the cropped canvas for this viewing angle.
            axis_3d.set_xlabel("")
            axis_3d.text2D(
                0.015,
                0.19,
                "X [in]",
                transform=axis_3d.transAxes,
                rotation=-20,
                ha="left",
                va="center",
            )
            axis_3d.set_ylabel("Y [in]", labelpad=18)
            axis_3d.set_zlabel("Z [in]", labelpad=5)
            axis_3d.xaxis.set_major_locator(MaxNLocator(nbins=4))
            axis_3d.yaxis.set_major_locator(MaxNLocator(nbins=7))
            axis_3d.zaxis.set_major_locator(MaxNLocator(nbins=5))
            axis_3d.tick_params(axis="both", which="major", labelsize=8, pad=1)
            axis_3d.set_title(
                f"Mode {mode} — {frequencies[mode]:.4f} Hz",
                fontsize=12.5,
                pad=4,
            )
            handles, labels_3d = axis_3d.get_legend_handles_labels()
            axis_3d.legend(
                handles,
                labels_3d,
                loc="lower center",
                bbox_to_anchor=(0.58, 0.015),
                ncol=2,
                fontsize=9,
                frameon=False,
            )

            # Tight limits with equal physical scale in inches on X, Y and Z.
            displayed_points = np.vstack((coordinate_array, deformed))
            lower = displayed_points.min(axis=0)
            upper = displayed_points.max(axis=0)
            padding = np.full(3, 10.0)
            plot_lower = lower - padding
            plot_upper = upper + padding
            plot_range = plot_upper - plot_lower
            axis_3d.set_xlim(plot_lower[0], plot_upper[0])
            axis_3d.set_ylim(plot_lower[1], plot_upper[1])
            axis_3d.set_zlim(plot_lower[2], plot_upper[2])
            axis_3d.set_box_aspect(tuple(plot_range), zoom=1.25)

            # Near-isometric orthographic view; positive X recedes into screen.
            axis_3d.set_proj_type("ortho")
            axis_3d.view_init(elev=25.0, azim=150.0)
            fig.subplots_adjust(left=0.08, right=0.97, bottom=0.03, top=0.89)
            fig.savefig(
                output / f"mode_{mode:02d}_mbdyn_3d.png",
                dpi=180,
                bbox_inches="tight",
                pad_inches=0.06,
            )
            plt.close(fig)
            print(f"Mode {mode}: MAC={mac:.9f}, RMS difference={rms:.3e}")

    with (output / "nastran_fem_vs_mbdyn_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "mode",
                "frequency_hz",
                "translational_MAC",
                "normalized_RMS_difference",
                "maximum_nodal_vector_difference",
            ]
        )
        writer.writerows(summary)
    print(f"Comparison diagrams written to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
