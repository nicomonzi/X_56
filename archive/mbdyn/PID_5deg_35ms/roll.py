"""Backward-compatible entry point for the roll dashboard."""

from plot_pid_results import main


if __name__ == "__main__":
    main(["--case", "roll"])
