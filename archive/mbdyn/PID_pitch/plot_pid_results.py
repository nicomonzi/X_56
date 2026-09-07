"""Generate English PID dashboards from the MBDyn NetCDF results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from cycler import cycler


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "output"
BASE_NODE = "990000"
PID_ID = "20"
TITLE_SIZE = 12
LABEL_SIZE = 14
LEGEND_SIZE = 12
TICK_SIZE = 11
THESIS_COLORS = [
    "#8B0000",  # dark red
    "#00008B",  # dark blue
    "#4B0082",  # dark purple
    "#006400",  # dark green
    "#CC5500",  # dark orange
]
TRACKING_ERROR_COLOR = "#87CEEB"  # light blue

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.weight": "normal",
        "axes.labelweight": "normal",
        "axes.titleweight": "normal",
        "axes.prop_cycle": cycler(color=THESIS_COLORS),
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)
CASES = {
    "roll": {
        "file": OUTPUT_DIR / "roll_10deg_40ms.nc",
        "axis": 0,
        "angle": "Roll",
        "rate": "Roll rate",
        "target": 10.0,
        "speed": 40.0,
        "rate_limit": None,
        "surface_limit": 30.0,
        "image": ROOT / "dashboard_roll_pid_10deg_40ms.png",
    },
    "pitch": {
        "file": OUTPUT_DIR / "pitch_5deg_30ms.nc",
        "axis": 1,
        "angle": "Pitch",
        "rate": "Pitch rate",
        "target": 5.0,
        "speed": 40.0,
        "rate_limit": None,
        "surface_limit": 10.0,
        "image": ROOT / "dashboard_pitch_pid_5deg_30ms.png",
    },
}


def component(dataset: nc.Dataset, variable: str, axis: int, steps: int) -> np.ndarray:
    """Read one component from either a 2-D or flattened MBDyn vector."""
    values = np.asarray(dataset.variables[variable][:], dtype=float)
    if values.ndim == 2:
        return values[:, axis]
    if values.size == steps * 3:
        return values[axis::3]
    if axis == 0 and values.size == steps:
        return values
    raise ValueError(f"Unexpected shape {values.shape} for '{variable}'")


def load_case(
    case_name: str, source: Path | None = None
) -> dict[str, np.ndarray]:
    config = CASES[case_name]
    input_file = source if source is not None else config["file"]
    with nc.Dataset(input_file, "r") as dataset:
        time = np.asarray(dataset.variables["time"][:], dtype=float)
        steps = time.size
        euler_name = f"node.struct.{BASE_NODE}.E"
        omega_name = f"node.struct.{BASE_NODE}.Omega"
        if euler_name not in dataset.variables:
            euler_name = f"node.struct.{BASE_NODE}.Phi"

        angle = np.degrees(component(dataset, euler_name, config["axis"], steps))
        rate = np.degrees(component(dataset, omega_name, config["axis"], steps))
        pid_output = np.degrees(
            np.asarray(dataset.variables[f"elem.loadable.{PID_ID}.output"][:], dtype=float)
        )
        rate_pid_name = "elem.loadable.21.output"
        if case_name == "pitch" and rate_pid_name in dataset.variables:
            pid_output += np.degrees(
                np.asarray(dataset.variables[rate_pid_name][:], dtype=float)
            )
        error = np.degrees(
            np.asarray(dataset.variables[f"elem.loadable.{PID_ID}.error"][:], dtype=float)
        )

    # The PID error is setpoint - measurement, so this reconstructs the exact
    # commanded angle from the MBDyn output without altering the data.
    target = angle + error
    return {
        "time": time,
        "angle": angle,
        "rate": rate,
        "pid_output": pid_output,
        "error": error,
        "target": target,
    }


def plot_case(
    case_name: str,
    source: Path | None = None,
    destination: Path | None = None,
) -> Path:
    config = CASES[case_name]
    data = load_case(case_name, source)
    time = data["time"]
    image = destination if destination is not None else config["image"]
    pdf = image.with_suffix(".pdf")

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(7.2, 10.5),
        sharex=True,
        constrained_layout=True,
    )
    fig.suptitle(
        rf"Target {config['angle']} angle = {config['target']:.0f}$^\circ$, "
        rf"$\mathbf{{u}}_\infty = \mathbf{{{config['speed']:.0f}\,m/s}}$",
        fontsize=TITLE_SIZE,
        fontweight="normal",
    )

    angle_symbol = r"$\phi$" if case_name == "roll" else r"$\theta$"
    rate_symbol = r"$p$" if case_name == "roll" else r"$q$"
    surface_symbol = r"$\delta_a$" if case_name == "roll" else r"$\delta_e$"

    axes[0].plot(
        time,
        data["angle"],
        color=THESIS_COLORS[0],
        linewidth=2.1,
        label=f"{config['angle']} Angle",
    )
    axes[0].plot(
        time,
        data["target"],
        color="0.15",
        linestyle="-",
        linewidth=1.2,
        label=f"Target {config['target']:.0f}°",
    )
    axes[0].fill_between(
        time,
        data["angle"],
        data["target"],
        color=TRACKING_ERROR_COLOR,
        alpha=0.45,
        label="Tracking Error",
    )
    axes[0].set_ylabel(f"{angle_symbol} [deg]", fontsize=LABEL_SIZE)

    axes[1].plot(
        time,
        data["rate"],
        color=THESIS_COLORS[0],
        linewidth=2.1,
        label=config["rate"],
    )
    axes[1].axhline(0.0, color="0.15", linestyle="-", linewidth=1.2)
    if config["rate_limit"] is not None:
        axes[1].axhline(
            config["rate_limit"],
            color="0.15",
            linestyle="-",
            linewidth=1.2,
            label=f"Rate constraint (±{config['rate_limit']:.0f}°/s)",
        )
        axes[1].axhline(
            -config["rate_limit"],
            color="0.15",
            linestyle="-",
            linewidth=1.2,
        )
    axes[1].set_ylabel(f"{rate_symbol} [deg/s]", fontsize=LABEL_SIZE)

    if case_name == "roll":
        axes[2].plot(
            time,
            -data["pid_output"],
            color=THESIS_COLORS[0],
            linewidth=1.8,
            label="Left Control Surfaces",
        )
        axes[2].plot(
            time,
            data["pid_output"],
            color=THESIS_COLORS[1],
            linewidth=1.8,
            label="Right Control Surfaces",
        )
    else:
        axes[2].plot(
            time,
            -data["pid_output"],
            color=THESIS_COLORS[0],
            linewidth=2.1,
            label="Control Surfaces",
        )
    axes[2].axhline(
        config["surface_limit"],
        color="0.15",
        linestyle="-",
        linewidth=1.2,
        label=f"Deflection limits ±{config['surface_limit']:.0f}°",
    )
    axes[2].axhline(
        -config["surface_limit"], color="0.15", linestyle="-", linewidth=1.2
    )
    axes[2].set_ylabel(f"{surface_symbol} [deg]", fontsize=LABEL_SIZE)

    for axis in axes:
        axis.set_xlim(time[0], time[-1])
        axis.set_xlabel(r"$t$ [s]", fontsize=LABEL_SIZE, fontweight="normal")
        axis.tick_params(
            axis="both",
            which="major",
            direction="in",
            labelsize=TICK_SIZE,
            length=5,
            width=0.8,
            top=True,
            right=True,
            labelbottom=True,
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
        axis.grid(
            True,
            which="major",
            linestyle="--",
            linewidth=0.6,
            alpha=0.35,
            color="0.45",
        )
        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
            spine.set_color("0.2")

        legend = axis.legend(
            loc="best",
            fontsize=LEGEND_SIZE,
            frameon=True,
            framealpha=0.90,
            edgecolor="0.75",
            fancybox=False,
        )
        for text in legend.get_texts():
            text.set_fontweight("normal")

    fig.savefig(image, dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    final_error = config["target"] - data["angle"][-1]
    print(
        f"{case_name.capitalize()}: final angle={data['angle'][-1]:.3f} deg, "
        f"final error={final_error:.3f} deg, "
        f"peak |deflection|={np.max(np.abs(data['pid_output'])):.3f} deg"
    )
    print(f"Saved {image}")
    print(f"Saved {pdf}")
    return image


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("roll", "pitch", "all"), default="all")
    parser.add_argument(
        "--input",
        type=Path,
        help="Read this NetCDF file (requires a single --case).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Save the dashboard to this path (requires a single --case).",
    )
    args = parser.parse_args(argv)
    if args.case == "all" and (args.input is not None or args.output is not None):
        parser.error("--input and --output require --case roll or --case pitch")
    selected = ("roll", "pitch") if args.case == "all" else (args.case,)
    for case_name in selected:
        plot_case(case_name, args.input, args.output)


if __name__ == "__main__":
    main()
