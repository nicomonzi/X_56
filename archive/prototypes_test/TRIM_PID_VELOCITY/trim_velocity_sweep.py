#!/usr/bin/env python3
"""Run an MBDyn longitudinal-trim sweep and create independent plots.

The script is intended to be placed in the same directory as the main MBDyn
input file and its relative ``include`` files.  For every airspeed it:

1. creates a temporary copy of the MBDyn input in the model directory;
2. replaces ``VELOCITY_MS`` without modifying the original file;
3. enables MBDyn text output in the temporary copy, so ``.jnt`` and ``.usr``
   can be read even when the original model uses ``no text``;
4. runs MBDyn with a separate output directory;
5. extracts Fz and My from total pin joint ``TRIM_JOINT`` and theta/delta_el
   from the saturated outputs of PID elements ``PID_THETA`` and
   ``PID_ELEVATOR``;
6. saves a CSV summary, a CSV history for every case, and six independent
   figures with English titles.

Default sweep: 20 m/s to 35 m/s, step 1 m/s.

Typical usage
-------------
    python3 trim_velocity_sweep.py --input main_trim.mbd

Optional executable path
------------------------
    python3 trim_velocity_sweep.py --input main_trim.mbd \
        --mbdyn /usr/local/mbdyn/bin/mbdyn

Dependencies
------------
    numpy, matplotlib
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


RAD_TO_DEG = 180.0 / math.pi
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167

DEFAULT_TRIM_JOINT = 23
DEFAULT_PID_THETA = 9101
DEFAULT_PID_ELEVATOR = 9102
DEFAULT_FINAL_TIME_S = 10.0
DEFAULT_TIME_STEP_S = 0.02
DEFAULT_SETTLING_WINDOW_S = 1.0


@dataclass
class CaseHistory:
    velocity_mps: float
    time_s: np.ndarray
    fz_n: np.ndarray
    my_nm: np.ndarray
    theta_deg: np.ndarray
    elevator_deg: np.ndarray


@dataclass
class CaseResult:
    velocity_mps: float
    theta_trim_deg: float
    elevator_trim_deg: float
    mean_fz_n: float
    mean_my_nm: float
    std_fz_n: float
    std_my_nm: float
    final_fz_n: float
    final_my_nm: float
    history: CaseHistory
    case_directory: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run an MBDyn trim sweep from 20 to 35 m/s and create plots."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("main_trim.mbd"),
        help="main MBDyn input file (default: main_trim.mbd)",
    )
    parser.add_argument(
        "--mbdyn",
        type=Path,
        default=None,
        help="path to the MBDyn executable; otherwise PATH/MBDYN_BIN is used",
    )
    parser.add_argument("--vmin", type=float, default=35.0, help="minimum airspeed [m/s]")
    parser.add_argument("--vmax", type=float, default=45.0, help="maximum airspeed [m/s]")
    parser.add_argument("--dv", type=float, default=1.0, help="airspeed increment [m/s]")
    parser.add_argument(
        "--settling-window",
        type=float,
        default=DEFAULT_SETTLING_WINDOW_S,
        help="final averaging window used for trim values [s]",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("trim_sweep_results"),
        help="result directory, relative to the model directory by default",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="continue with the next airspeed if one MBDyn case fails",
    )
    return parser.parse_args()


def locate_mbdyn(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit.expanduser())

    env_value = os.environ.get("MBDYN_BIN")
    if env_value:
        candidates.append(Path(env_value).expanduser())

    found = shutil.which("mbdyn")
    if found:
        candidates.append(Path(found))

    candidates.extend(
        [
            Path("/usr/local/mbdyn/bin/mbdyn"),
            Path("/usr/local/bin/mbdyn"),
            Path("/usr/bin/mbdyn"),
        ]
    )

    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()

    raise FileNotFoundError(
        "MBDyn executable not found. Use --mbdyn /path/to/mbdyn or set MBDYN_BIN."
    )


def fortran_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def build_velocity_grid(vmin: float, vmax: float, dv: float) -> np.ndarray:
    if vmin <= 0.0:
        raise ValueError("vmin must be positive")
    if vmax < vmin:
        raise ValueError("vmax must be greater than or equal to vmin")
    if dv <= 0.0:
        raise ValueError("dv must be positive")

    count = int(math.floor((vmax - vmin) / dv + 1.0e-10)) + 1
    values = vmin + dv * np.arange(count, dtype=float)

    if values[-1] < vmax - 1.0e-9:
        values = np.append(values, vmax)
    else:
        values[-1] = vmax

    return values


def read_integer_constant(source: str, name: str, fallback: int) -> int:
    pattern = re.compile(
        rf"set\s*:\s*(?:ifndef\s+)?const\s+integer\s+{re.escape(name)}\s*=\s*([+-]?\d+)\s*;",
        flags=re.IGNORECASE,
    )
    match = pattern.search(source)
    return fallback if match is None else int(match.group(1))


def read_real_card(source: str, card_name: str, fallback: float) -> float:
    pattern = re.compile(
        rf"\b{re.escape(card_name)}\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)\s*;",
        flags=re.IGNORECASE,
    )
    match = pattern.search(source)
    return fallback if match is None else fortran_float(match.group(1))


def replace_velocity(source: str, velocity_mps: float) -> str:
    pattern = re.compile(
        r"(set\s*:\s*(?:ifndef\s+)?const\s+real\s+VELOCITY_MS\s*=\s*)([^;]+)(;)",
        flags=re.IGNORECASE,
    )
    replaced, count = pattern.subn(
        lambda match: f"{match.group(1)}{velocity_mps:.16g}{match.group(3)}",
        source,
        count=1,
    )
    if count != 1:
        raise RuntimeError(
            "Exactly one declaration 'set: const real VELOCITY_MS = ...;' is required."
        )
    return replaced


def enable_text_output(source: str) -> str:
    """Enable ASCII output in a generated copy without touching the original."""

    pattern = re.compile(
        r"(output\s+results\s*:\s*netcdf\b)([^;]*)(;)",
        flags=re.IGNORECASE | re.DOTALL,
    )

    def replacement(match: re.Match[str]) -> str:
        options = match.group(2)
        options = re.sub(
            r",\s*no\s+text\b", ", text", options, flags=re.IGNORECASE
        )
        if not re.search(r"(?:^|,)\s*text\b", options, flags=re.IGNORECASE):
            options = f"{options}, text"
        return f"{match.group(1)}{options}{match.group(3)}"

    return pattern.sub(replacement, source, count=1)


def prepare_case_input(
    base_source: str,
    model_directory: Path,
    case_directory: Path,
    velocity_mps: float,
) -> Path:
    generated = enable_text_output(replace_velocity(base_source, velocity_mps))
    name = f"__trim_sweep_V{velocity_mps:06.2f}.mbd"
    case_input = model_directory / name
    case_input.write_text(generated, encoding="utf-8")
    shutil.copy2(case_input, case_directory / "case_input.mbd")
    return case_input


def times_from_lines(lines: Iterable[str]) -> np.ndarray:
    values: list[float] = []
    active = False

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        header = line.lstrip("# ").split()
        if "Step" in header and "Time" in header and "TStep" in header:
            active = True
            continue

        if not active or line.startswith("#"):
            continue

        fields = line.replace(",", " ").split()
        if len(fields) < 2:
            continue

        try:
            int(fields[0])
            time_value = fortran_float(fields[1])
        except ValueError:
            continue

        if math.isfinite(time_value):
            values.append(time_value)

    return np.asarray(values, dtype=float)


def infer_times(sample_count: int, final_time_s: float, time_step_s: float) -> np.ndarray:
    if sample_count < 1:
        return np.asarray([], dtype=float)

    expected_intervals = int(round(final_time_s / time_step_s))
    if sample_count == expected_intervals + 1:
        return np.linspace(0.0, final_time_s, sample_count)
    if sample_count == expected_intervals:
        return np.arange(1, sample_count + 1, dtype=float) * time_step_s

    return np.arange(sample_count, dtype=float) * time_step_s


def read_joint_reactions(path: Path, joint_label: int) -> tuple[np.ndarray, np.ndarray]:
    """Read local Fz and My from a total pin joint .jnt file."""

    fz_lbf: list[float] = []
    my_lbfin: list[float] = []

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            fields = raw.split()
            if len(fields) < 7:
                continue

            try:
                values = [fortran_float(field) for field in fields]
            except ValueError:
                continue

            if int(round(values[0])) != joint_label:
                continue

            # Total pin joint text output:
            # label, local Fx, Fy, Fz, local Mx, My, Mz, ...
            fz_lbf.append(values[3])
            my_lbfin.append(values[5])

    if not fz_lbf:
        raise RuntimeError(f"joint {joint_label} not found in {path}")

    return np.asarray(fz_lbf, dtype=float), np.asarray(my_lbfin, dtype=float)


def read_pid_outputs(
    path: Path,
    theta_pid_label: int,
    elevator_pid_label: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Read saturated PID outputs from column 3 of the .usr file."""

    theta_rad: list[float] = []
    elevator_rad: list[float] = []

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            fields = raw.split()
            if len(fields) < 3:
                continue

            try:
                label = int(round(fortran_float(fields[0])))
                saturated_output = fortran_float(fields[2])
            except ValueError:
                continue

            if label == theta_pid_label:
                theta_rad.append(saturated_output)
            elif label == elevator_pid_label:
                elevator_rad.append(saturated_output)

    if not theta_rad:
        raise RuntimeError(f"PID_THETA {theta_pid_label} not found in {path}")
    if not elevator_rad:
        raise RuntimeError(f"PID_ELEVATOR {elevator_pid_label} not found in {path}")

    return np.asarray(theta_rad, dtype=float), np.asarray(elevator_rad, dtype=float)


def align_histories(
    time_s: np.ndarray,
    fz_lbf: np.ndarray,
    my_lbfin: np.ndarray,
    theta_rad: np.ndarray,
    elevator_rad: np.ndarray,
    velocity_mps: float,
) -> CaseHistory:
    sample_count = min(
        len(time_s), len(fz_lbf), len(my_lbfin), len(theta_rad), len(elevator_rad)
    )
    if sample_count < 3:
        raise RuntimeError("not enough output samples to build time histories")

    return CaseHistory(
        velocity_mps=velocity_mps,
        time_s=time_s[-sample_count:].copy(),
        fz_n=fz_lbf[-sample_count:].copy() * LBF_TO_N,
        my_nm=my_lbfin[-sample_count:].copy() * LBFIN_TO_NM,
        theta_deg=theta_rad[-sample_count:].copy() * RAD_TO_DEG,
        elevator_deg=elevator_rad[-sample_count:].copy() * RAD_TO_DEG,
    )


def write_history_csv(path: Path, history: CaseHistory) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "time_s",
                "Fz_N",
                "My_Nm",
                "theta_deg",
                "delta_elevator_deg",
            ]
        )
        writer.writerows(
            zip(
                history.time_s,
                history.fz_n,
                history.my_nm,
                history.theta_deg,
                history.elevator_deg,
            )
        )


def run_case(
    mbdyn: Path,
    model_directory: Path,
    base_source: str,
    output_root: Path,
    velocity_mps: float,
    trim_joint: int,
    pid_theta: int,
    pid_elevator: int,
    final_time_s: float,
    time_step_s: float,
    settling_window_s: float,
) -> CaseResult:
    case_directory = output_root / f"V_{velocity_mps:05.1f}_mps"
    if case_directory.exists():
        shutil.rmtree(case_directory)
    case_directory.mkdir(parents=True)

    case_input = prepare_case_input(
        base_source, model_directory, case_directory, velocity_mps
    )
    output_prefix = case_directory / "result"

    command = [
        str(mbdyn),
        "-f",
        case_input.name,
        "-o",
        str(output_prefix.resolve()),
    ]

    try:
        completed = subprocess.run(
            command,
            cwd=model_directory,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    finally:
        case_input.unlink(missing_ok=True)

    console_path = case_directory / "console.log"
    console_path.write_text(completed.stdout, encoding="utf-8", errors="replace")

    if completed.returncode != 0:
        tail = "\n".join(completed.stdout.splitlines()[-40:])
        raise RuntimeError(
            f"MBDyn failed at V={velocity_mps:g} m/s with return code "
            f"{completed.returncode}.\nLast console lines:\n{tail}"
        )

    out_path = output_prefix.with_suffix(".out")
    jnt_path = output_prefix.with_suffix(".jnt")
    usr_path = output_prefix.with_suffix(".usr")

    missing = [str(path) for path in (jnt_path, usr_path) if not path.exists()]
    if missing:
        available = ", ".join(sorted(path.name for path in case_directory.iterdir()))
        raise RuntimeError(
            "Required MBDyn text output is missing: "
            + ", ".join(missing)
            + f". Available files: {available}"
        )

    fz_lbf, my_lbfin = read_joint_reactions(jnt_path, trim_joint)
    theta_rad, elevator_rad = read_pid_outputs(usr_path, pid_theta, pid_elevator)

    time_s = np.asarray([], dtype=float)
    if out_path.exists():
        time_s = times_from_lines(
            out_path.read_text(encoding="utf-8", errors="replace").splitlines()
        )
    if time_s.size == 0:
        time_s = times_from_lines(completed.stdout.splitlines())
    if time_s.size == 0:
        sample_count = min(len(fz_lbf), len(my_lbfin), len(theta_rad), len(elevator_rad))
        time_s = infer_times(sample_count, final_time_s, time_step_s)

    history = align_histories(
        time_s,
        fz_lbf,
        my_lbfin,
        theta_rad,
        elevator_rad,
        velocity_mps,
    )
    write_history_csv(case_directory / "history.csv", history)

    mask = history.time_s >= history.time_s[-1] - settling_window_s - 1.0e-12
    if np.count_nonzero(mask) < 3:
        mask = np.arange(len(history.time_s)) >= max(0, len(history.time_s) - 10)

    return CaseResult(
        velocity_mps=velocity_mps,
        theta_trim_deg=float(np.mean(history.theta_deg[mask])),
        elevator_trim_deg=float(np.mean(history.elevator_deg[mask])),
        mean_fz_n=float(np.mean(history.fz_n[mask])),
        mean_my_nm=float(np.mean(history.my_nm[mask])),
        std_fz_n=float(np.std(history.fz_n[mask])),
        std_my_nm=float(np.std(history.my_nm[mask])),
        final_fz_n=float(history.fz_n[-1]),
        final_my_nm=float(history.my_nm[-1]),
        history=history,
        case_directory=case_directory,
    )


def write_summary_csv(output_root: Path, results: Sequence[CaseResult]) -> Path:
    path = output_root / "trim_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "velocity_mps",
                "theta_trim_deg",
                "delta_elevator_trim_deg",
                "mean_Fz_N",
                "mean_My_Nm",
                "std_Fz_N",
                "std_My_Nm",
                "final_Fz_N",
                "final_My_Nm",
                "case_directory",
            ]
        )
        for result in sorted(results, key=lambda item: item.velocity_mps):
            writer.writerow(
                [
                    f"{result.velocity_mps:.12g}",
                    f"{result.theta_trim_deg:.12g}",
                    f"{result.elevator_trim_deg:.12g}",
                    f"{result.mean_fz_n:.12g}",
                    f"{result.mean_my_nm:.12g}",
                    f"{result.std_fz_n:.12g}",
                    f"{result.std_my_nm:.12g}",
                    f"{result.final_fz_n:.12g}",
                    f"{result.final_my_nm:.12g}",
                    str(result.case_directory),
                ]
            )
    return path


def finish_axes(ax: plt.Axes, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)


def save_velocity_plot(
    output_root: Path,
    velocities: np.ndarray,
    values: np.ndarray,
    filename: str,
    title: str,
    ylabel: str,
) -> Path:
    path = output_root / filename
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(velocities, values, marker="o", linewidth=1.6)
    finish_axes(ax, title, "Airspeed [m/s]", ylabel)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def save_history_plot(
    output_root: Path,
    results: Sequence[CaseResult],
    selector: str,
    filename: str,
    title: str,
    ylabel: str,
    zero_line: bool,
) -> Path:
    path = output_root / filename
    fig, ax = plt.subplots(figsize=(9.2, 5.6))

    for result in sorted(results, key=lambda item: item.velocity_mps):
        values = getattr(result.history, selector)
        ax.plot(
            result.history.time_s,
            values,
            linewidth=1.15,
            label=f"{result.velocity_mps:g} m/s",
        )

    if zero_line:
        ax.axhline(0.0, linewidth=0.9)

    finish_axes(ax, title, "Time [s]", ylabel)
    ax.legend(
        title="Airspeed",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def create_plots(output_root: Path, results: Sequence[CaseResult]) -> list[Path]:
    ordered = sorted(results, key=lambda item: item.velocity_mps)
    velocities = np.asarray([result.velocity_mps for result in ordered])
    theta_trim = np.asarray([result.theta_trim_deg for result in ordered])
    elevator_trim = np.asarray([result.elevator_trim_deg for result in ordered])

    paths = [
        save_velocity_plot(
            output_root,
            velocities,
            theta_trim,
            "theta_trim_vs_velocity.png",
            "Trim Pitch Angle versus Airspeed",
            r"Pitch angle $\theta$ [deg]",
        ),
        save_velocity_plot(
            output_root,
            velocities,
            elevator_trim,
            "elevator_trim_vs_velocity.png",
            "Trim Elevator Deflection versus Airspeed",
            r"Elevator deflection $\delta_{el}$ [deg]",
        ),
        save_history_plot(
            output_root,
            ordered,
            "fz_n",
            "vertical_reaction_time_histories.png",
            "Vertical Reaction Force Convergence",
            r"Vertical reaction $F_z$ [N]",
            True,
        ),
        save_history_plot(
            output_root,
            ordered,
            "my_nm",
            "pitching_moment_time_histories.png",
            "Pitching Moment Convergence",
            r"Pitching moment $M_y$ [N m]",
            True,
        ),
        save_history_plot(
            output_root,
            ordered,
            "theta_deg",
            "pitch_angle_time_histories.png",
            "Pitch Angle Time Histories",
            r"Pitch angle $\theta$ [deg]",
            False,
        ),
        save_history_plot(
            output_root,
            ordered,
            "elevator_deg",
            "elevator_deflection_time_histories.png",
            "Elevator Deflection Time Histories",
            r"Elevator deflection $\delta_{el}$ [deg]",
            False,
        ),
    ]
    return paths


def main() -> int:
    args = parse_args()

    input_path = args.input.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"MBDyn input file not found: {input_path}")
    if args.settling_window <= 0.0:
        raise ValueError("settling-window must be positive")

    model_directory = input_path.parent
    output_root = (
        args.output_dir.expanduser()
        if args.output_dir.is_absolute()
        else model_directory / args.output_dir
    ).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    base_source = input_path.read_text(encoding="utf-8", errors="replace")
    velocities = build_velocity_grid(args.vmin, args.vmax, args.dv)
    mbdyn = locate_mbdyn(args.mbdyn)

    trim_joint = read_integer_constant(base_source, "TRIM_JOINT", DEFAULT_TRIM_JOINT)
    pid_theta = read_integer_constant(base_source, "PID_THETA", DEFAULT_PID_THETA)
    pid_elevator = read_integer_constant(
        base_source, "PID_ELEVATOR", DEFAULT_PID_ELEVATOR
    )
    final_time_s = read_real_card(base_source, "final time", DEFAULT_FINAL_TIME_S)
    time_step_s = read_real_card(base_source, "time step", DEFAULT_TIME_STEP_S)

    print(f"MBDyn executable: {mbdyn}")
    print(f"Input file: {input_path}")
    print(f"Output directory: {output_root}")
    print(
        "Airspeeds [m/s]: " + ", ".join(f"{velocity:g}" for velocity in velocities)
    )
    print(
        f"Labels: TRIM_JOINT={trim_joint}, PID_THETA={pid_theta}, "
        f"PID_ELEVATOR={pid_elevator}"
    )

    results: list[CaseResult] = []
    failures: list[tuple[float, str]] = []

    for index, velocity_mps in enumerate(velocities, start=1):
        print(f"\n[{index}/{len(velocities)}] Running V = {velocity_mps:g} m/s")
        try:
            result = run_case(
                mbdyn=mbdyn,
                model_directory=model_directory,
                base_source=base_source,
                output_root=output_root,
                velocity_mps=float(velocity_mps),
                trim_joint=trim_joint,
                pid_theta=pid_theta,
                pid_elevator=pid_elevator,
                final_time_s=final_time_s,
                time_step_s=time_step_s,
                settling_window_s=args.settling_window,
            )
        except (OSError, RuntimeError, ValueError) as exc:
            message = str(exc)
            failures.append((float(velocity_mps), message))
            print(f"FAILED: {message}", file=sys.stderr)
            if not args.continue_on_error:
                raise
            continue

        results.append(result)
        print(
            f"theta={result.theta_trim_deg:+.6f} deg, "
            f"delta_el={result.elevator_trim_deg:+.6f} deg, "
            f"mean Fz={result.mean_fz_n:+.6e} N, "
            f"mean My={result.mean_my_nm:+.6e} N m"
        )

    if not results:
        raise RuntimeError("No MBDyn case completed successfully")

    summary_path = write_summary_csv(output_root, results)
    plot_paths = create_plots(output_root, results)

    if failures:
        failure_path = output_root / "failures.csv"
        with failure_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["velocity_mps", "error"])
            writer.writerows(failures)
        print(f"Failure report: {failure_path}")

    print("\nCreated files:")
    print(f"  {summary_path}")
    for path in plot_paths:
        print(f"  {path}")

    return 2 if failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)