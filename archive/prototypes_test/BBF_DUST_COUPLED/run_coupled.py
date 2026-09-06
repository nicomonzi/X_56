#!/usr/bin/env python3
"""Launch one coupled X-56 MBDyn/DUST case; no post-processing is performed."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent
DT = 0.02
MBDYN_PYTHON_PATHS = (
    Path("/home/nicomonzi/src/mbdyn/libraries/libmbc"),
    Path("/home/nicomonzi/src/mbdyn/libraries/libmbc/.libs"),
    Path("/usr/local/mbdyn/libexec/mbpy"),
)


def executable(value: str | None, default: str) -> Path:
    candidate = value or shutil.which(default)
    if not candidate:
        raise RuntimeError(f"Required executable not found: {default}")
    path = Path(candidate).resolve()
    if not path.is_file():
        raise RuntimeError(f"Executable does not exist: {path}")
    return path


def dust_has_precice(dust_bin: Path) -> bool:
    """Reject the known non-coupled build before starting any solver."""
    result = subprocess.run(
        ["strings", str(dust_bin)],
        check=True,
        capture_output=True,
        text=True,
    )
    return "compiled without #USE_PRECICE" not in result.stdout


def make_case_inputs(velocity: float, case_dir: Path) -> tuple[Path, Path]:
    mbd_text = (ROOT / "main_bbf_dust.mbd").read_text()
    mbd_text = mbd_text.replace(
        "set: const real V_INF = 30.0;",
        f"set: const real V_INF = {velocity:.12g};",
        1,
    )
    # The generated case lives below output/. Resolve includes explicitly so
    # MBDyn never searches relative to the result directory (also on /mnt/c).
    mbd_text = mbd_text.replace(
        "../BBF_AUTO_TRIM/INCLUDE/",
        f"{(ROOT.parent / 'BBF_AUTO_TRIM' / 'INCLUDE').resolve()}/",
    )
    mbd_text = mbd_text.replace(
        "./INCLUDE/",
        f"{(ROOT / 'INCLUDE').resolve()}/",
    )
    mbd_input = case_dir / "case_input.mbd"
    mbd_input.write_text(mbd_text)

    speed_ips = velocity / 0.0254
    dust_text = (ROOT / "dust" / "dust.in").read_text()
    dust_text = dust_text.replace(
        "u_inf = (/1181.1023622047244, 0.0, 0.0/)",
        f"u_inf = (/{speed_ips:.15g}, 0.0, 0.0/)",
        1,
    )
    dust_output = (case_dir / "dust" / "case").resolve()
    dust_text = dust_text.replace(
        "basename = ./Output/case",
        f"basename = {dust_output}",
        1,
    )
    dust_text = dust_text.replace(
        "geometry_file = geo_input.h5",
        f"geometry_file = {(ROOT / 'dust' / 'geo_input.h5').resolve()}",
        1,
    )
    dust_text = dust_text.replace(
        "reference_file = References.in",
        f"reference_file = {(ROOT / 'dust' / 'References.in').resolve()}",
        1,
    )
    dust_input = case_dir / "dust_input.in"
    dust_input.write_text(dust_text)
    return mbd_input, dust_input


def stop_process(process: subprocess.Popen | None) -> None:
    if process is not None and process.poll() is None:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


def protect_or_archive_previous_result(case_dir: Path) -> None:
    """Refuse completed data; preserve and move aside an interrupted run."""
    result = case_dir / "case.nc"
    if not result.exists():
        return

    complete = False
    try:
        from netCDF4 import Dataset

        with Dataset(result) as dataset:
            time_values = dataset.variables["time"]
            complete = len(time_values) > 0 and float(time_values[-1]) >= 29.74
    except Exception:
        # A damaged/unfinished NetCDF file is archived like any interrupted run.
        complete = False

    if complete:
        raise SystemExit(f"Complete result already exists: {result}")

    archive = case_dir / f"incomplete_attempt_{int(time.time())}"
    archive.mkdir(parents=True)
    for path in case_dir.glob("case.*"):
        shutil.move(str(path), archive / path.name)
    dust_archive = archive / "dust"
    for path in (case_dir / "dust").glob("case*"):
        dust_archive.mkdir(exist_ok=True)
        shutil.move(str(path), dust_archive / path.name)
    print(f"Archived interrupted result in: {archive}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the coupled X-56 MBDyn/DUST aeroelastic model."
    )
    parser.add_argument("--velocity", type=float, default=30.0, help="m/s")
    parser.add_argument(
        "--dust-bin",
        default=str(ROOT / "build_dust" / "bin" / "dust"),
        help="preCICE-enabled DUST executable",
    )
    parser.add_argument(
        "--dust-pre-bin",
        default=str(ROOT / "build_dust" / "bin" / "dust_pre"),
        help="preCICE-enabled DUST preprocessor executable",
    )
    parser.add_argument("--mbdyn-bin", default="/usr/local/mbdyn/bin/mbdyn")
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="validate files and executables without starting the solvers",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.velocity <= 0.0:
        raise SystemExit("Velocity must be positive")

    dust_bin = executable(args.dust_bin, "dust")
    dust_pre_bin = executable(args.dust_pre_bin, "dust_pre")
    mbdyn_bin = executable(args.mbdyn_bin, "mbdyn")
    if not dust_has_precice(dust_bin):
        raise SystemExit(
            f"{dust_bin} was compiled without preCICE support.\n"
            "Rebuild DUST with WITH_PRECICE=YES and pass that binary with "
            "--dust-bin. No uncoupled run was started."
        )

    case_dir = ROOT / "output" / f"V_{args.velocity:06.2f}_mps"
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "dust").mkdir(exist_ok=True)
    if not args.check_only:
        protect_or_archive_previous_result(case_dir)
    mbd_input, dust_input = make_case_inputs(args.velocity, case_dir)

    # Build the DUST geometry once from the original parametric X-56 mesh.
    subprocess.run(
        [str(dust_pre_bin), "dust_pre.in"],
        cwd=ROOT / "dust",
        check=True,
    )
    if args.check_only:
        print(f"Coupled inputs validated for {args.velocity:g} m/s")
        return

    # preCICE leaves this handshake directory behind after an interrupted or
    # failed start. Reusing it makes the next run abort before data exchange.
    # It contains connection metadata only, never aerodynamic/structural data.
    precice_state = ROOT / "precice-run"
    if precice_state.exists():
        shutil.rmtree(precice_state)

    for path in MBDYN_PYTHON_PATHS:
        if path.is_dir():
            sys.path.insert(0, str(path))
    from mbdyn_interface import MBDynInterface
    from coupling_adapter import CouplingAdapter

    dust_log = (case_dir / "dust.log").open("w")
    mbdyn_log = (case_dir / "mbdyn.log").open("w")
    dust_process = None
    mbdyn_process = None
    with tempfile.TemporaryDirectory(prefix="mbdyn_dust_") as temp_dir:
        socket_path = str(Path(temp_dir) / "mbdyn.sock")
        environment = os.environ.copy()
        environment["MBSOCK"] = socket_path
        try:
            dust_process = subprocess.Popen(
                [str(dust_bin), str(dust_input)],
                cwd=ROOT,
                env=environment,
                stdout=dust_log,
                stderr=subprocess.STDOUT,
            )
            mbdyn_process = subprocess.Popen(
                [
                    str(mbdyn_bin),
                    "-f", str(mbd_input),
                    "-o", str(case_dir / "case"),
                ],
                cwd=ROOT,
                env=environment,
                stdout=mbdyn_log,
                stderr=subprocess.STDOUT,
            )

            # MBDyn creates the Unix socket; wait briefly before negotiating.
            deadline = time.monotonic() + 30.0
            while not Path(socket_path).exists():
                if mbdyn_process.poll() is not None:
                    raise RuntimeError(
                        f"MBDyn stopped early; inspect {case_dir / 'mbdyn.log'}"
                    )
                if time.monotonic() > deadline:
                    raise RuntimeError("Timed out waiting for the MBDyn socket")
                time.sleep(0.05)

            mbdyn = MBDynInterface(
                socket_path, str(ROOT / "structure" / "refConfigNodes.in")
            )
            CouplingAdapter(
                mbdyn, str(ROOT / "precice-config.xml")
            ).run(DT)

            mbdyn_code = mbdyn_process.wait(timeout=30)
            dust_code = dust_process.wait(timeout=30)
            if mbdyn_code or dust_code:
                raise RuntimeError(
                    f"Solver exit codes: MBDyn={mbdyn_code}, DUST={dust_code}"
                )
        finally:
            stop_process(mbdyn_process)
            stop_process(dust_process)
            mbdyn_log.close()
            dust_log.close()

    print(f"Coupled result: {case_dir / 'case.nc'}")


if __name__ == "__main__":
    main()
