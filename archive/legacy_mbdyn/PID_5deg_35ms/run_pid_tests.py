"""Run the two 35 m/s MBDyn PID cases and generate their dashboards."""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
MBDYN = "/usr/local/mbdyn/bin/mbdyn"
CASES = (
    ("main_roll.mbd", "output/roll_10deg_40ms"),
    ("main_pitch.mbd", "output/pitch_5deg_30ms"),
)


def main() -> None:
    (ROOT / "output").mkdir(exist_ok=True)
    for input_file, output_prefix in CASES:
        print(f"Running {input_file}...")
        subprocess.run(
            [MBDYN, "-f", input_file, "-o", output_prefix],
            cwd=ROOT,
            check=True,
        )
    subprocess.run(
        [sys.executable, "plot_pid_results.py", "--case", "all"],
        cwd=ROOT,
        check=True,
    )


if __name__ == "__main__":
    main()
