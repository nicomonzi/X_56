#!/usr/bin/env python3
"""Portable CHECK/SMOKE/PRODUCTION launcher for the X-56A DUST case."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile
import time

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent
VELOCITY = 60.8421
SETTINGS = {
    "smoke": {"mesh": "COARSE", "dt": 0.002, "final": 0.10, "steps": 50},
    "production": {"mesh": "FINE", "dt": 0.002, "final": 9.50, "steps": 4750},
}
PORTABLE_SCAN = ("run_case.py", "model", "meshes", "adapters", "tools",
                 "config/machine.env.example", "README_RUN.md", "reference")


def load_machine_env() -> None:
    path = ROOT / "config/machine.env"
    if not path.exists():
        return
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        name, value = line.split("=", 1)
        os.environ.setdefault(name.strip(), value.strip())


def resolve_executable(variable: str, command: str) -> Path:
    value = os.environ.get(variable) or shutil.which(command)
    if not value:
        raise RuntimeError(f"{variable} is unset and {command!r} is not on PATH")
    path = Path(value).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise RuntimeError(f"Invalid {variable}: {path}")
    return path


def detect_precice_major() -> int:
    requested = os.environ.get("PRECICE_MAJOR")
    if requested:
        return int(requested)
    script = "import precice; print(3 if hasattr(precice, 'Participant') else 2)"
    result = subprocess.run([sys.executable, "-c", script], capture_output=True,
                            text=True, timeout=15, env=os.environ.copy())
    if result.returncode:
        raise RuntimeError("Cannot import preCICE Python bindings:\n" + result.stderr.strip())
    return int(result.stdout.strip().splitlines()[-1])


def set_python_paths() -> None:
    value = os.environ.get("MBDYN_PYTHON_PATH", "")
    paths = [str(Path(item).expanduser()) for item in value.split(os.pathsep) if item]
    for path in reversed(paths):
        if path not in sys.path:
            sys.path.insert(0, path)
    if paths:
        old = os.environ.get("PYTHONPATH", "")
        os.environ["PYTHONPATH"] = os.pathsep.join(paths + ([old] if old else []))


def solver_environment(threads: int) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update({
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": str(threads),
        "OMP_PLACES": "cores",
        "OMP_PROC_BIND": "close",
    })
    return environment


def print_resources(threads: int, environment: dict[str, str]) -> None:
    print(f"hostname: {socket.gethostname()}")
    print("date: " + time.strftime("%Y-%m-%d %H:%M:%S %Z"))
    try:
        print("uptime: " + Path("/proc/uptime").read_text().split()[0] + " s")
        mem = next(line for line in Path("/proc/meminfo").read_text().splitlines()
                   if line.startswith("MemAvailable:"))
        print("available RAM: " + " ".join(mem.split()[1:]))
        load1 = os.getloadavg()[0]
        print(f"load average (1 min): {load1:.2f}")
        if load1 > max(threads, (os.cpu_count() or threads) * 0.75):
            print("WARNING: machine is already heavily loaded")
    except (OSError, StopIteration):
        print("system load/RAM: unavailable")
    print(f"requested threads (hard maximum): {threads}")
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_PLACES", "OMP_PROC_BIND"):
        print(f"{name}={environment[name]}")


def write_dust_pre(level: str, output: Path) -> Path:
    path = output / "dust_pre.in"
    path.write_text(
        "comp_name = X56\n"
        f"geo_file = meshes/{level}/parametric_mesh.in\n"
        "ref_tag = centerbody\n"
        f"file_name = {output.relative_to(ROOT)}/geo_input.h5\n"
    )
    return path


def generate_geometry(level: str, output: Path, dust_pre: Path) -> dict:
    config = write_dust_pre(level, output)
    subprocess.run([str(dust_pre), str(config.relative_to(ROOT))], cwd=ROOT, check=True)
    geometry_file = output / "geo_input.h5"
    with h5py.File(geometry_file) as h5:
        geometry = h5["Components/Comp001/Geometry"]
        result = {"elements": len(geometry["ee"]), "points": len(geometry["rr"]),
                  "coupling_nodes": len(geometry["CouplingNodes"])}
    return result


def make_inputs(kind: str, output: Path, precice_major: int) -> tuple[Path, Path, Path]:
    main_template = ROOT / f"model/mbdyn/main_{kind}.mbd"
    text = main_template.read_text().replace(
        "./INCLUDE/", f"{(ROOT/'model/mbdyn/INCLUDE').resolve()}/")
    mbdyn_input = output / f"{kind}.mbd"
    mbdyn_input.write_text(text)

    config_source = ROOT / ("model/precice-config-v3.xml" if precice_major >= 3
                            else "model/precice-config.xml")
    config_text = config_source.read_text()
    if kind == "smoke":
        config_text = config_text.replace('max-time value="9.50"', 'max-time value="0.10"')
    precice_config = output / "precice-config.xml"
    precice_config.write_text(config_text)

    dust_text = (ROOT / f"model/dust/dust_{kind}.in").read_text()
    prefix = f"output/{kind}"
    dust_text = dust_text.replace(f"{prefix}/geo_input.h5", str((output/"geo_input.h5").resolve()))
    dust_text = dust_text.replace(f"{prefix}/precice-config.xml", str(precice_config.resolve()))
    dust_text = dust_text.replace(f"{prefix}/dust/case", str((output/"dust/case").resolve()))
    dust_text = dust_text.replace("model/dust/References.in", str((ROOT/"model/dust/References.in").resolve()))
    dust_input = output / f"dust_{kind}.in"
    dust_input.write_text(dust_text)
    return mbdyn_input, dust_input, precice_config


def print_modes() -> None:
    print("Retained provisional modal basis:")
    with (ROOT / "reference/modal_basis.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            print(f"  FEM {int(row['fem_mode']):2d}: {float(row['frequency_hz']):8.5f} Hz  {row['description']}")


def portable_path_audit() -> list[str]:
    hits = []
    for entry in PORTABLE_SCAN:
        start = ROOT / entry
        paths = [start] if start.is_file() else start.rglob("*")
        for path in paths:
            if not path.is_file() or path.suffix in {".fem", ".h5", ".vtu", ".png"}:
                continue
            try: text = path.read_text()
            except UnicodeDecodeError: continue
            forbidden = ("/" + "home/", "/" + "mnt/c/")
            if any(marker in text for marker in forbidden):
                hits.append(str(path.relative_to(ROOT)))
    return hits


def production_blockers() -> list[str]:
    path = ROOT / "reports/fine_mesh_audit.json"
    if not path.exists():
        return ["FINE mesh audit report is missing"]
    return json.loads(path.read_text()).get("blocking_reasons", [])


def probe_mbdyn_parse(binary: Path, input_file: Path, output: Path,
                      environment: dict[str, str]) -> None:
    """Parse through the external-force card, then stop before any time step."""
    socket_path = output / "check_mbdyn.sock"
    environment = environment.copy(); environment["MBSOCK"] = str(socket_path)
    process = subprocess.Popen([str(binary), "-f", str(input_file), "-o", str(output/"parse")],
                               cwd=ROOT, env=environment, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True)
    deadline = time.monotonic() + 8
    while process.poll() is None and not socket_path.exists() and time.monotonic() < deadline:
        time.sleep(.05)
    stop(process)
    stdout = process.stdout.read() if process.stdout else ""
    (output/"mbdyn_parse.log").write_text(stdout)
    if not socket_path.exists() and "Reading Force(9500)" not in stdout:
        raise RuntimeError("MBDyn input did not parse through external force; see mbdyn_parse.log")


def probe_dust_parse(binary: Path, input_file: Path, output: Path,
                     environment: dict[str, str]) -> None:
    """Let DUST reach its preCICE wait, then stop before a coupled time step."""
    process = subprocess.Popen([str(binary), str(input_file)], cwd=output,
                               env=environment, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True)
    deadline = time.monotonic() + 5
    while process.poll() is None and time.monotonic() < deadline:
        time.sleep(.05)
    was_alive = process.poll() is None
    stop(process)
    stdout = process.stdout.read() if process.stdout else ""
    (output/"dust_parse.log").write_text(stdout)
    parsed = "Reading input parameters from file" in stdout and "ERROR in \"dust main\"" not in stdout
    if not (was_alive or parsed):
        raise RuntimeError("DUST input parse probe failed; see dust_parse.log")


def check_only(args, bins, precice_major: int, environment: dict[str, str]) -> None:
    output = ROOT / "output/check"
    (output / "dust").mkdir(parents=True, exist_ok=True)
    mesh = generate_geometry("COARSE", output, bins["dust_pre"])
    mbd, dust, config = make_inputs("smoke", output, precice_major)
    required = [mbd, dust, config, ROOT/"model/mbdyn/INCLUDE/mbdyn_modal.fem",
                ROOT/"model/dust/coupling_nodes.in", ROOT/"model/structure/refConfigNodes.in"]
    missing = [str(path) for path in required if not path.is_file()]
    if missing: raise RuntimeError("Missing files: " + ", ".join(missing))
    if not os.access(output, os.W_OK): raise RuntimeError(f"Output is not writable: {output}")
    if importlib.util.find_spec("numpy") is None or importlib.util.find_spec("h5py") is None:
        raise RuntimeError("numpy/h5py are unavailable")
    if importlib.util.find_spec("mbc_py_interface") is None:
        raise RuntimeError("mbc_py_interface is unavailable; set MBDYN_PYTHON_PATH")
    if portable_path_audit():
        raise RuntimeError("Machine-specific paths remain in portable files: " +
                           ", ".join(portable_path_audit()))
    if "compiled without #USE_PRECICE" in subprocess.run(
            ["strings", str(bins["dust"])], capture_output=True, text=True, check=True).stdout:
        raise RuntimeError("DUST executable was compiled without preCICE")
    probe_mbdyn_parse(bins["mbdyn"], mbd, output, environment)
    probe_dust_parse(bins["dust"], dust, output, environment)
    production_output = ROOT / "output/check_production"
    (production_output / "dust").mkdir(parents=True, exist_ok=True)
    production_mesh = generate_geometry("FINE", production_output, bins["dust_pre"])
    production_mbd, production_dust, _ = make_inputs(
        "production", production_output, precice_major)
    probe_mbdyn_parse(bins["mbdyn"], production_mbd, production_output, environment)
    probe_dust_parse(bins["dust"], production_dust, production_output, environment)
    print(f"CHECK PASS: executables, Python requirements, paths, writable output, and {mesh}")
    print(f"FINE production input parse geometry: {production_mesh}")
    print(f"preCICE Python major API detected: {precice_major}")
    print("MBDyn and DUST parsed to their coupling waits; no physical step advanced.")
    print_modes()
    blockers = production_blockers()
    print("PRODUCTION READY: NO")
    for blocker in blockers: print("  BLOCKER: " + blocker)


def stop(process: subprocess.Popen | None) -> None:
    if process is not None and process.poll() is None:
        process.terminate()
        try: process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill(); process.wait()


def postprocess(kind: str, output: Path, dust_post: Path) -> None:
    # output_start=F: DUST stores dt_out, 2*dt_out, ... strictly before tend.
    frames = 4 if kind == "smoke" else 474
    directory = output / "paraview"; directory.mkdir(exist_ok=True)
    config = output / "dust_post.in"
    config.write_text(f"""data_basename = {(output/'dust/case').resolve()}
basename = {(directory/'x56').resolve()}
analysis = {{
  type = viz
  name = coupled
  start_res = 1
  end_res = {frames}
  step_res = 1
  format = vtk
  wake = T
  separate_wake = T
  variable = velocity
  variable = vorticity
  variable = pressure
  component = all
}}
""")
    subprocess.run([str(dust_post), str(config)], cwd=ROOT, check=True,
                   stdout=(output/"dust_post.log").open("w"), stderr=subprocess.STDOUT)
    subprocess.run([sys.executable, str(ROOT/"tools/fix_vtu_xml.py"), str(directory)], check=True)
    for suffix in ("", "_wpan", "_wpart"):
        lines = ['<?xml version="1.0"?>',
                 '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
                 '  <Collection>']
        for frame in range(1, frames+1):
            lines.append(f'    <DataSet timestep="{frame*0.02:.5f}" group="" part="0" file="x56_coupled{suffix}-{frame:04d}.vtu"/>')
        lines += ['  </Collection>', '</VTKFile>']
        (directory/f"x56_coupled{suffix}.pvd").write_text("\n".join(lines)+"\n")


def coupled_run(kind: str, args, bins, precice_major: int,
                environment: dict[str, str]) -> None:
    if kind == "production":
        blockers = production_blockers()
        if blockers:
            raise SystemExit("PRODUCTION BLOCKED:\n- " + "\n- ".join(blockers))
    settings = SETTINGS[kind]
    output = ROOT / f"output/{kind}"
    if output.exists() and any(output.glob("dust/case_res_*.h5")):
        raise SystemExit(f"Existing result protected: {output}")
    (output / "dust").mkdir(parents=True, exist_ok=True)
    mesh = generate_geometry(settings["mesh"], output, bins["dust_pre"])
    mbd_input, dust_input, config = make_inputs(kind, output, precice_major)
    print(f"{kind.upper()} mesh/settings: {mesh}; dt={settings['dt']}; "
          f"final={settings['final']}; steps={settings['steps']}")
    if kind == "smoke": print("SMOKE TEST - NOT PHYSICALLY VALID")
    print_resources(args.threads, environment)

    sys.path.insert(0, str(ROOT / "adapters"))
    from mbdyn_interface import MBDynInterface
    from coupling_adapter import CouplingAdapter
    dust_process = mbdyn_process = None
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix=f"x56_{kind}_") as temp:
        socket_path = str(Path(temp) / "mbdyn.sock")
        environment["MBSOCK"] = socket_path
        with (output/"dust.log").open("w") as dust_log, (output/"mbdyn.log").open("w") as mbdyn_log:
            try:
                dust_process = subprocess.Popen([str(bins["dust"]), str(dust_input)],
                    cwd=ROOT, env=environment, stdout=dust_log, stderr=subprocess.STDOUT)
                mbdyn_process = subprocess.Popen([str(bins["mbdyn"]), "-f", str(mbd_input),
                    "-o", str(output/kind)], cwd=ROOT, env=environment,
                    stdout=mbdyn_log, stderr=subprocess.STDOUT)
                deadline = time.monotonic() + 60
                while not Path(socket_path).exists():
                    if mbdyn_process.poll() is not None:
                        raise RuntimeError(f"MBDyn stopped early; inspect {output/'mbdyn.log'}")
                    if time.monotonic() > deadline: raise RuntimeError("MBDyn socket timeout")
                    time.sleep(.05)
                interface = MBDynInterface(socket_path, str(ROOT/"model/structure/refConfigNodes.in"))
                def state_at(t):
                    if kind == "smoke": return "SMOKE"
                    if t < 5.5: return "CAPTURE_TRIM"
                    if t < 6.0: return "READY"
                    if t < 6.85: return "RAP"
                    if t < 9.40: return "OPEN_LOOP"
                    if t < 9.50: return "RECOVERY"
                    return "END"
                CouplingAdapter(interface, str(config),
                                str(output/"coupled_response.csv"),
                                state_function=state_at).run(settings["dt"])
                codes = (mbdyn_process.wait(timeout=60), dust_process.wait(timeout=60))
                if codes != (0, 0): raise RuntimeError(f"Solver exit codes: {codes}")
            finally:
                stop(mbdyn_process); stop(dust_process)
    summary = {"classification": "SMOKE TEST - NOT PHYSICALLY VALID" if kind == "smoke" else "production",
               "wall_time_s": time.monotonic()-start, "threads": args.threads, **settings, **mesh}
    (output/"run_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    subprocess.run([
        sys.executable, str(ROOT/"tools/extract_mbdyn_text.py"), str(output),
        "--prefix", kind,
    ], check=True)
    postprocess(kind, output, bins["dust_post"])
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check-only", action="store_true")
    group.add_argument("--smoke", action="store_true")
    group.add_argument("--production", action="store_true")
    parser.add_argument("--threads", type=int, choices=(8, 12, 16), default=12)
    args = parser.parse_args()
    load_machine_env(); set_python_paths()
    bins = {"mbdyn": resolve_executable("MBDYN_BIN", "mbdyn"),
            "dust": resolve_executable("DUST_BIN", "dust"),
            "dust_pre": resolve_executable("DUST_PRE_BIN", "dust_pre"),
            "dust_post": resolve_executable("DUST_POST_BIN", "dust_post")}
    major = detect_precice_major()
    environment = solver_environment(args.threads)
    print_resources(args.threads, environment)
    for name, path in bins.items(): print(f"{name.upper()}={path}")
    if args.check_only: check_only(args, bins, major, environment)
    elif args.smoke: coupled_run("smoke", args, bins, major, environment)
    else: coupled_run("production", args, bins, major, environment)


if __name__ == "__main__":
    main()
