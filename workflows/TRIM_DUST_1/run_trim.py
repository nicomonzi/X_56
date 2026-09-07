#!/usr/bin/env python3
"""Coupled DUST/MBDyn trim of the X-56A at 55 m/s (no SAS, no BFF)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent
CASE = ROOT / "case"
WORK = ROOT / "work" / "current"
RESULTS = ROOT / "results"
DT = 0.004
TRIM_TIME = 8.0
VELOCITY_MPS = 55.0
VINF_INPS = VELOCITY_MPS / 0.0254
WEIGHT_LBF = 419.4399625
CG = np.array([163.187383385809, 0.110529571088, 101.239797358848])
MESH_ORIGIN = np.array([100.0, 0.0, 100.16000370])
CANONICAL_SHA256 = "c448d6ffec0fdde60aa7459b3ae631bf5f9954c44709810655845caed39df7b3"
WRITE_FIELDS = ("Position", "Velocity", "Rotation", "AngularVelocity")
READ_FIELDS = ("Force", "Moment")


def executable(env_name: str, preferred: Path, fallback: str) -> Path:
    candidates = []
    if os.environ.get(env_name):
        candidates.append(Path(os.environ[env_name]))
    candidates.append(preferred)
    found = shutil.which(fallback)
    if found:
        candidates.append(Path(found))
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise RuntimeError(f"eseguibile mancante: imposta {env_name}")


def binaries() -> dict[str, Path]:
    patched = Path.home() / "dust-patched" / "build-user" / "bin"
    dust = executable("DUST_BIN", patched / "dust", "dust")
    dust_pre = executable("DUST_PRE_BIN", dust.with_name("dust_pre"), "dust_pre")
    if dust.parent != dust_pre.parent:
        raise RuntimeError("dust e dust_pre devono provenire dalla stessa build accoppiata")
    return {
        "dust": dust,
        "dust_pre": dust_pre,
        "mbdyn": executable("MBDYN_BIN", Path("/usr/local/mbdyn/bin/mbdyn"), "mbdyn"),
    }


def render(text: str, values: dict[str, object]) -> str:
    for key, value in values.items():
        text = text.replace(f"__{key}__", str(value))
    unresolved = sorted(set(re.findall(r"__([A-Z0-9_]+)__", text)))
    if unresolved:
        raise RuntimeError(f"placeholder non risolti: {unresolved}")
    return text


def render_mesh() -> str:
    source = (CASE / "mesh_canonical.in").read_bytes()
    digest = hashlib.sha256(source).hexdigest()
    if digest != CANONICAL_SHA256:
        raise RuntimeError(
            "mesh_canonical.in non coincide con la mesh DUST_MESH validata; "
            f"SHA256={digest}")
    lines = source.decode().splitlines()
    output: list[str] = []
    hinge_index = 0
    inserted = False
    skipping_function = False
    chord_syntax = os.environ.get("DUST_CHORD_SYNTAX", "cosine_le")
    if chord_syntax not in ("cosine_le", "cosineLE"):
        raise RuntimeError("DUST_CHORD_SYNTAX deve essere cosine_le oppure cosineLE")
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("hinge_tag"):
            hinge_index += 1
        if skipping_function:
            if stripped == "}":
                skipping_function = False
            continue
        if stripped == "el_type = p" and not inserted:
            output.append(line)
            output.extend((
                "", "! Coupling acts only on the 39 elastic wing stations.",
                "coupled = T", "coupling_type = rbf",
                "coupling_node_file = coupling_nodes.in",
            ))
            inserted = True
            continue
        if stripped.startswith("type_chord"):
            line = f"type_chord = {chord_syntax}"
        elif stripped.startswith("hinge_rotation_input"):
            first = 40 + 2 * (hinge_index - 1)
            output.extend((
                "hinge_rotation_input = coupling",
                "hinge_rotation_coupling = {",
                "coupling_node_subset = range",
                f"coupling_node_first = {first}",
                f"coupling_node_last = {first + 1}",
                "}",
            ))
            continue
        elif stripped.startswith("hinge_rotation_function"):
            skipping_function = True
            continue
        output.append(line)
    if not inserted:
        raise RuntimeError("header della mesh canonica non riconosciuto")
    text = "\n".join(output) + "\n"
    if text.count("hinge_rotation_input = coupling") != 10:
        raise RuntimeError("la mesh runtime non contiene esattamente 10 hinge accoppiate")
    return text


def prepare(duration: float = TRIM_TIME) -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    (WORK / "dust").mkdir(parents=True)
    os.symlink(CASE / "airfoilsection", WORK / "airfoilsection", target_is_directory=True)
    for name in ("References.in", "dust_pre.in"):
        shutil.copy2(CASE / name, WORK / name)
    structural_nodes = np.loadtxt(CASE / "coupling_nodes_59.in").reshape(-1, 3)
    np.savetxt(WORK / "coupling_nodes.in", structural_nodes - MESH_ORIGIN,
               fmt="%.12f")
    (WORK / "mesh.in").write_text(render_mesh())
    (WORK / "precice.xml").write_text(render(
        (CASE / "precice_trim.xml").read_text(),
        {"FINAL_TIME": f"{duration:.9g}"}))
    dust = render((CASE / "dust.in.in").read_text(), {
        "PRECICE_XML": "precice.xml", "FINAL_TIME": f"{duration:.9g}"})
    (WORK / "dust.in").write_text(dust)
    mbdyn = render((CASE / "trim.mbd.in").read_text(), {
        "FINAL_TIME": f"{duration:.9g}",
        "CASE_DIR": str(CASE)})
    (WORK / "trim.mbd").write_text(mbdyn)


def validate_geometry(path: Path) -> tuple[int, int, int, int]:
    with h5py.File(path) as h5:
        group = h5["Components/Comp001"]
        geometry = group["Geometry"]
        shape = (geometry["ee"].shape[0], geometry["rr"].shape[0],
                 geometry["CouplingNodes"].shape[0], int(group["Hinges/n_hinges"][()]))
    expected = (2440, 2583, 59, 10)
    if shape != expected:
        raise RuntimeError(f"mesh inattesa {shape}; attesa {expected}")
    return shape


class MBDynSocket:
    def __init__(self, path: str):
        mbpy = Path(os.environ.get("MBDYN_MBPY", "/usr/local/mbdyn/libexec/mbpy"))
        if mbpy.is_dir() and str(mbpy) not in sys.path:
            sys.path.insert(0, str(mbpy))
        from mbc_py_interface import mbcNodal

        self.reference = np.loadtxt(CASE / "coupling_nodes_59.in").reshape(-1, 3)
        self.precice_reference = self.reference - MESH_ORIGIN
        count = len(self.reference)
        self.nodal = mbcNodal(path, "", 0, -1, 0, 1, 0, count, 0, 0x100, 1)
        self.nodal.negotiate()
        self.nodal.recv()
        self.data = {
            "Position": self._array(self.nodal.n_x),
            "Velocity": self._array(self.nodal.n_xp),
            "Rotation": self._array(self.nodal.n_theta),
            "AngularVelocity": self._array(self.nodal.n_omega),
            "Force": np.zeros((count, 3)), "Moment": np.zeros((count, 3)),
        }

    def _array(self, values):
        return np.asarray(values, dtype=float).reshape(len(self.reference), 3).copy()

    def receive(self) -> bool:
        if self.nodal.recv():
            return True
        self.data["Position"] = self._array(self.nodal.n_x)
        self.data["Velocity"] = self._array(self.nodal.n_xp)
        self.data["Rotation"] = self._array(self.nodal.n_theta)
        self.data["AngularVelocity"] = self._array(self.nodal.n_omega)
        return False

    def send(self, converged: bool, scale: float) -> bool:
        self.nodal.n_f[:] = (scale * self.data["Force"]).reshape(-1)
        self.nodal.n_m[:] = (scale * self.data["Moment"]).reshape(-1)
        return bool(self.nodal.send(converged))

    def close(self):
        self.nodal.destroy()


def run_adapter(socket_path: str, duration: float) -> None:
    import precice

    mbdyn = MBDynSocket(socket_path)
    participant = precice.Participant("MBDyn", str(WORK / "precice.xml"), 0, 1)
    mesh_name = "MBDynNodes"
    ids = participant.set_mesh_vertices(mesh_name, mbdyn.precice_reference)

    def write_data():
        for field in WRITE_FIELDS:
            invalid = np.argwhere(~np.isfinite(mbdyn.data[field]))
            if invalid.size:
                node_index, component = invalid[0]
                raise RuntimeError(
                    f"MBDyn ha prodotto NaN/Inf nel campo {field}, "
                    f"nodo coupling {int(node_index) + 1}, componente {int(component) + 1}")
            values = (mbdyn.data[field] - MESH_ORIGIN
                      if field == "Position" else mbdyn.data[field])
            participant.write_data(mesh_name, field, ids, values)

    with (WORK / "coupling.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time_s", "left_tip_dz_in", "right_tip_dz_in",
                         "symmetric_tip_dz_in", "total_fx_lbf", "total_fy_lbf",
                         "total_fz_lbf", "global_mx_cg_lbfin", "global_my_cg_lbfin",
                         "global_mz_cg_lbfin"))
        if participant.requires_initial_data():
            write_data()
        participant.initialize()
        print("preCICE connesso: trim accoppiato in avanzamento", flush=True)
        current_time = 0.0
        try:
            while participant.is_coupling_ongoing():
                ramp_x = min(max(current_time / 0.40, 0.0), 1.0)
                load_scale = ramp_x * ramp_x * (3.0 - 2.0 * ramp_x)
                if mbdyn.send(False, load_scale) or mbdyn.receive():
                    break
                write_data()
                step = min(DT, participant.get_max_time_step_size())
                participant.advance(step)
                if not participant.is_coupling_ongoing():
                    break
                for field in READ_FIELDS:
                    mbdyn.data[field] = participant.read_data(mesh_name, field, ids, step)
                if mbdyn.send(True, load_scale) or mbdyn.receive():
                    break
                current_time += step
                displacement = mbdyn.data["Position"] - mbdyn.reference
                force = mbdyn.data["Force"].sum(axis=0)
                moment = np.sum(mbdyn.data["Moment"] + np.cross(
                    mbdyn.data["Position"] - CG, mbdyn.data["Force"]), axis=0)
                writer.writerow((f"{current_time:.9f}", displacement[19, 2],
                                 displacement[38, 2],
                                 0.5 * (displacement[19, 2] + displacement[38, 2]),
                                 *force, *moment))
                stream.flush()
                if not np.isfinite(displacement).all() or not np.isfinite(force).all():
                    raise RuntimeError(f"stato non finito a t={current_time:.3f} s")
        finally:
            participant.finalize()
            mbdyn.close()


def terminate(process) -> None:
    if process is not None and process.poll() is None:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)


def terminal_complete(duration: float) -> bool:
    try:
        with (WORK / "coupling.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        return (float(rows[-1]["time_s"]) >= duration - DT - 1.e-9
                and "Computations Finished" in (WORK / "dust.log").read_text())
    except (FileNotFoundError, IndexError, KeyError, ValueError):
        return False


def coupled_run(duration: float, threads: int) -> None:
    bins = binaries()
    prepare(duration)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS=str(threads),
               OMP_PLACES="cores", OMP_PROC_BIND="close")
    with (WORK / "dust_pre.log").open("w") as log:
        completed = subprocess.run([str(bins["dust_pre"]), "dust_pre.in"], cwd=WORK,
                                   env=env, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError("dust_pre fallito: controlla work/current/dust_pre.log")
    validate_geometry(WORK / "geo_input.h5")
    socket_path = f"/tmp/zeno_trim55_{os.getpid()}.sock"
    Path(socket_path).unlink(missing_ok=True)
    env["MBSOCK"] = socket_path
    dust_process = mbdyn_process = None
    with (WORK / "dust.log").open("w") as dust_log, (WORK / "mbdyn.log").open("w") as mbdyn_log:
        try:
            dust_process = subprocess.Popen([str(bins["dust"]), "dust.in"], cwd=WORK,
                                            env=env, stdout=dust_log,
                                            stderr=subprocess.STDOUT, start_new_session=True)
            mbdyn_process = subprocess.Popen(
                [str(bins["mbdyn"]), "-f", "trim.mbd", "-o", str(WORK / "case")],
                cwd=WORK, env=env, stdout=mbdyn_log, stderr=subprocess.STDOUT,
                start_new_session=True)
            for _ in range(900):
                if Path(socket_path).exists():
                    break
                if dust_process.poll() is not None or mbdyn_process.poll() is not None:
                    raise RuntimeError("solver fermo prima del coupling: controlla i log")
                time.sleep(0.1)
            else:
                raise RuntimeError("timeout in attesa del socket MBDyn")
            previous = Path.cwd()
            try:
                os.chdir(WORK)
                run_adapter(socket_path, duration)
            finally:
                os.chdir(previous)
            codes = (mbdyn_process.wait(60), dust_process.wait(60))
            if codes != (0, 0) and not (codes[1] == 0 and terminal_complete(duration)):
                raise RuntimeError(f"codici solver MBDyn={codes[0]} DUST={codes[1]}")
        finally:
            terminate(mbdyn_process)
            terminate(dust_process)
            Path(socket_path).unlink(missing_ok=True)


def trim(args) -> None:
    """Run one online trim, freeze its commands, then verify for two seconds."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    if args.analyse_existing:
        print("analisi della run esistente: nessun solver viene rilanciato")
    else:
        print("single coupled run: trim 0.4--6.0 s, frozen verification 6.0--8.0 s")
        coupled_run(TRIM_TIME, args.threads)
    with (WORK / "coupling.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    frozen = [row for row in rows if float(row["time_s"]) >= 6.0]
    if not frozen:
        raise RuntimeError("finestra congelata 6--8 s assente")
    times = np.array([float(row["time_s"]) for row in frozen])
    fz_values = np.array([float(row["total_fz_lbf"]) for row in frozen])
    my_values = np.array([float(row["global_my_cg_lbfin"]) for row in frozen])
    tip_values = np.array([float(row["symmetric_tip_dz_in"]) for row in frozen])
    fz = float(np.mean(fz_values)); my = float(np.mean(my_values))
    # MBDyn writes classic NetCDF. scipy is already part of the numerical
    # Python stack used by the case and avoids requiring the optional netCDF4
    # module on the server.
    from scipy.io import netcdf_file
    with netcdf_file(WORK / "case.nc", "r", mmap=False) as dataset:
        mb_time = np.asarray(dataset.variables["time"].data, dtype=float)
        mask = mb_time >= 6.0
        if not np.any(mask):
            raise RuntimeError("finestra MBDyn congelata 6--8 s assente")
        reaction_fz = np.asarray(dataset.variables["elem.joint.23.F"].data[:, 2])[mask]
        reaction_my = np.asarray(dataset.variables["elem.joint.23.M"].data[:, 1])[mask]
        pitch = np.asarray(dataset.variables["node.struct.990000.Phi"].data[:, 1])[mask]
        bfl = np.asarray(dataset.variables["elem.joint.1004.Phi"].data[:, 1])[mask]
        bfr = np.asarray(dataset.variables["elem.joint.2004.Phi"].data[:, 1])[mask]
        modal = np.asarray(dataset.variables["elem.joint.5.a"].data)[mask]
    result = {
        "velocity_mps": VELOCITY_MPS,
        "method": "single-run online DUST/MBDyn trim with coupled native hinges",
        "dt_s": DT, "duration_s": TRIM_TIME,
        "active_trim_s": [0.4, 6.0], "frozen_verification_s": [6.0, 8.0],
        "mesh_panels": 2440, "coupling_nodes": 59,
        "trim_controls": "pitch and symmetric BFL/BFR; WF1-WF4 held at zero",
        "alpha_deg": float(np.rad2deg(np.mean(pitch))),
        "bfl_deg": float(np.rad2deg(np.mean(bfl))),
        "bfr_deg": float(np.rad2deg(np.mean(bfr))),
        "mean_fz_lbf": fz, "target_fz_lbf": WEIGHT_LBF,
        "mean_my_lbfin": my,
        "mean_trim_joint_fz_lbf": float(np.mean(reaction_fz)),
        "mean_trim_joint_my_lbfin": float(np.mean(reaction_my)),
        "trim_joint_fz_slope_lbf_s": float(np.polyfit(mb_time[mask], reaction_fz, 1)[0]),
        "trim_joint_my_slope_lbfin_s": float(np.polyfit(mb_time[mask], reaction_my, 1)[0]),
        "fz_slope_lbf_s": float(np.polyfit(times, fz_values, 1)[0]),
        "my_slope_lbfin_s": float(np.polyfit(times, my_values, 1)[0]),
        "tip_slope_in_s": float(np.polyfit(times, tip_values, 1)[0]),
        "modal_q": {str(mode): float(np.mean(modal[:, column]))
                    for column, mode in enumerate(range(7, 19))},
    }
    result["accepted"] = bool(
        abs(result["mean_trim_joint_fz_lbf"]) <= 5.0
        and abs(result["mean_trim_joint_my_lbfin"]) <= 30.0
        and abs(result["trim_joint_fz_slope_lbf_s"]) <= 5.0
        and abs(result["trim_joint_my_slope_lbfin_s"]) <= 30.0)
    result_path = RESULTS / "trim_55.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    print(f"result: {result_path} accepted={result['accepted']}")
    if not result["accepted"] and not args.analyse_existing:
        raise RuntimeError("trim non assestato nella finestra congelata 6--8 s")


def check(args) -> None:
    bins = binaries()
    prepare(TRIM_TIME)
    env = os.environ.copy(); env["OMP_NUM_THREADS"] = str(args.threads)
    with (WORK / "dust_pre.log").open("w") as log:
        code = subprocess.run([str(bins["dust_pre"]), "dust_pre.in"], cwd=WORK,
                              env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    if code:
        raise RuntimeError("dust_pre fallito: controlla work/current/dust_pre.log")
    shape = validate_geometry(WORK / "geo_input.h5")
    print(f"CHECK OK: panels={shape[0]} points={shape[1]} coupling_nodes={shape[2]} hinges={shape[3]}")
    print(f"DUST={bins['dust']}\nMBDyn={bins['mbdyn']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="preprocessa e valida soltanto la mesh")
    parser.add_argument("--analyse-existing", action="store_true",
                        help="analizza work/current senza rilanciare DUST o MBDyn")
    parser.add_argument("--threads", type=int, default=24)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("threads deve essere positivo")
    (ROOT / "work").mkdir(exist_ok=True); RESULTS.mkdir(exist_ok=True)
    check(args) if args.check else trim(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAILED: {error}", file=sys.stderr)
        raise SystemExit(1)
