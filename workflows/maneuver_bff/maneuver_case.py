#!/usr/bin/env python3
"""Render and run one time-domain maneuver X-56 BFF identification case."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CONFIG = json.loads((ROOT / "study_config.json").read_text())
BASELINE_VALUE = Path(CONFIG["baseline_directory"])
BASELINE = BASELINE_VALUE if BASELINE_VALUE.is_absolute() else (ROOT / BASELINE_VALUE).resolve()
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")
MODEL_REVISION = "pullup_surface_hold_v3"
ANGLE_MODEL_REVISION = "pullup_angle_surface_hold_v2_new_altitude"
DIVE_PULLUP_MODEL_REVISION = "dive_pullup_surface_hold_v4_n1_dt01"


def model_revision(point: "ManeuverPoint") -> str:
    if point.family == "pullup_angle":
        return ANGLE_MODEL_REVISION
    if point.family == "dive_pullup":
        return DIVE_PULLUP_MODEL_REVISION
    return MODEL_REVISION


def _load_baseline_module():
    # run_case imports helper modules that live beside it.
    if str(BASELINE) not in sys.path:
        sys.path.insert(0, str(BASELINE))
    spec = importlib.util.spec_from_file_location("bff_open_loop_run_case", BASELINE / "run_case.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("impossibile caricare BFF_open_loop/run_case.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_BASE = None


def baseline_module():
    global _BASE
    if _BASE is None:
        _BASE = _load_baseline_module()
    return _BASE


@dataclass(frozen=True)
class ManeuverPoint:
    family: str
    velocity_mps: float
    load_factor: float
    bank_angle_deg: float = 0.0
    excited: bool = True
    pitch_angle_deg: float = 0.0
    pitch_rate_deg_s: float = 0.0
    nominal_load_factor: float = 1.0

    def validate(self) -> None:
        if self.family not in {"pullup", "pullup_angle", "dive_pullup", "roll"}:
            raise ValueError("family deve essere pullup, pullup_angle, dive_pullup o roll")
        if not 50.0 <= self.velocity_mps <= 70.0:
            raise ValueError("la ROM e' validata solo tra 50 e 70 m/s TAS")
        if not 1.0 <= self.load_factor <= 2.5:
            raise ValueError("load_factor deve essere compreso tra 1 e 2.5")
        if self.family in {"pullup", "pullup_angle", "dive_pullup"} and abs(self.bank_angle_deg) > 1e-12:
            raise ValueError("il pull-up frozen deve essere simmetrico, bank_angle=0")
        if self.family not in {"pullup_angle", "dive_pullup"} and abs(self.pitch_angle_deg) > 1e-12:
            raise ValueError("pitch_angle_deg e' ammesso soltanto per pullup_angle/dive_pullup")
        if self.family == "pullup_angle" and not 0.0 <= self.pitch_angle_deg <= 12.0:
            raise ValueError("pullup_angle richiede 0 <= pitch_angle <= 12 deg")
        if self.family == "dive_pullup" and not 0.0 <= self.pitch_angle_deg <= 15.0:
            raise ValueError("dive_pullup richiede 0 <= pitch_angle <= 15 deg")
        if self.family != "dive_pullup" and abs(self.pitch_rate_deg_s) > 1e-12:
            raise ValueError("pitch_rate_deg_s e' ammesso soltanto per dive_pullup")
        if self.family == "dive_pullup" and not 0.0 <= self.pitch_rate_deg_s <= 10.0:
            raise ValueError("dive_pullup richiede 0 <= pitch_rate <= 10 deg/s")
        if not 1.0 <= self.nominal_load_factor <= 2.5:
            raise ValueError("nominal_load_factor deve essere compreso tra 1 e 2.5")
        if abs(self.bank_angle_deg) > 45.0:
            raise ValueError("il primo studio roll e' limitato a |bank| <= 45 deg")

    @property
    def stem(self) -> str:
        pitch_tag = (
            f"_pitch_{self.pitch_angle_deg:+06.2f}"
            if self.family in {"pullup_angle", "dive_pullup"} else ""
        )
        rate_tag = (
            f"_qcmd_{self.pitch_rate_deg_s:+06.2f}_nnom_{self.nominal_load_factor:05.2f}"
            if self.family == "dive_pullup" else ""
        )
        value = (
            f"{self.family}_V_{self.velocity_mps:07.3f}_n_{self.load_factor:05.3f}"
            f"_bank_{self.bank_angle_deg:+06.2f}{pitch_tag}{rate_tag}"
            f"_{'excited' if self.excited else 'shadow'}"
        )
        return value.replace(".", "p").replace("+", "p").replace("-", "m")


def _replace_real(text: str, name: str, value: float) -> str:
    text, count = re.subn(
        rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*[^;]+;",
        f"set: const real {name} = {value:.10e};",
        text,
    )
    if count != 1:
        raise RuntimeError(f"{name}: attesa una definizione, trovate {count}")
    return text


def _replace_pid_setpoint(text: str, pid_name: str, replacement: str) -> str:
    pattern = rf"(user defined:\s*{re.escape(pid_name)},\s*pid,.*?setpoint,)\s*const,\s*[^,]+,"
    text, count = re.subn(pattern, rf"\1 {replacement}", text, count=1, flags=re.S)
    if count != 1:
        raise RuntimeError(f"setpoint {pid_name} non trovato")
    return text


def resolved_timing(point: ManeuverPoint) -> dict[str, float]:
    """Return timing, extending aggressive dive--pull-ups automatically."""
    if point.family == "pullup_angle":
        return dict(CONFIG["angle_timing"])
    if point.family != "dive_pullup":
        return dict(CONFIG["timing"])
    timing = dict(CONFIG["dive_pullup_timing"])
    if point.pitch_rate_deg_s <= 1e-12:
        # The zero-command control case is steady longitudinal flight.  Keep
        # the 0.01 s step validated by BFF_open_loop: at 0.02 s the rigid-body
        # controller develops a numerical oscillation before SAS release,
        # so a case merely labelled n=1 would not actually remain at 1 g.
        timing["time_step_s"] = 0.01
        return timing
    amplitude = point.pitch_angle_deg
    release_lead = float(timing["pullup_release_lead_s"])
    sas_duration = float(timing["sas_off_duration_s"])
    dive_ramp = max(
        float(timing["minimum_dive_ramp_s"]),
        amplitude / float(timing["maximum_dive_rate_deg_s"]),
    )
    pullup_ramp = 2.0 * amplitude / point.pitch_rate_deg_s
    expected_ramp = release_lead + sas_duration
    if abs(pullup_ramp - expected_ramp) > 1e-6:
        raise ValueError("ampiezza e pitch-rate non centrano la finestra SAS-off")
    return_ramp = max(
        float(timing["minimum_return_ramp_s"]),
        amplitude / float(timing["maximum_return_rate_deg_s"]),
    )
    start = float(timing["maneuver_start_s"])
    dive_end = start + dive_ramp
    pullup_start = dive_end + float(timing["dive_hold_s"])
    pullup_end = pullup_start + pullup_ramp
    maneuver_end = pullup_end + float(timing["recovery_hold_s"]) + return_ramp
    control_end = maneuver_end + float(timing["control_recovery_s"])
    timing.update({
        "maneuver_rise_s": dive_ramp,
        "dive_ramp_s": dive_ramp,
        "pullup_ramp_s": pullup_ramp,
        "return_ramp_s": return_ramp,
        "sas_off_start_s": pullup_start + release_lead,
        "final_time_s": control_end + float(timing["final_settling_s"]),
    })
    return timing


def render_point(point: ManeuverPoint) -> str:
    point.validate()
    timing = resolved_timing(point)
    text = baseline_module().render_case(point.velocity_mps)
    text = _replace_real(text, "FINAL_TIME", float(timing["final_time_s"]))
    if "time_step_s" in timing:
        text = _replace_real(text, "TIME_STEP", float(timing["time_step_s"]))
    text = _replace_real(text, "SETTLING_END", float(timing["sas_off_start_s"]) - 0.5)
    text = _replace_real(text, "SAS_OFF_START", float(timing["sas_off_start_s"]))
    text = _replace_real(text, "SAS_OFF_DURATION", float(timing["sas_off_duration_s"]))

    maneuver_start = float(timing["maneuver_start_s"])
    maneuver_rise = float(timing["maneuver_rise_s"])
    maneuver_return = float(timing.get("maneuver_return_start_s", timing["final_time_s"] + 1.0))
    maneuver_fall = float(timing.get("maneuver_fall_s", maneuver_rise))
    angle_enable = float(point.family == "pullup_angle" and point.pitch_angle_deg > 1e-12)
    dive_pullup_enable = float(point.family == "dive_pullup" and point.pitch_angle_deg > 1e-12)
    dive_end = maneuver_start + float(timing.get("dive_ramp_s", maneuver_rise))
    pullup_start = dive_end + float(timing.get("dive_hold_s", 0.0))
    pullup_end = pullup_start + float(timing.get("pullup_ramp_s", maneuver_rise))
    recovery_end = pullup_end + float(timing.get("recovery_hold_s", 0.0))
    dive_maneuver_end = recovery_end + float(timing.get("return_ramp_s", maneuver_fall))
    dive_control_recovery_end = dive_maneuver_end + float(timing.get("control_recovery_s", 2.0))
    maneuver_enable = float(
        point.load_factor > 1.0 + 1e-12
        or abs(point.bank_angle_deg) > 1e-12
        or angle_enable > 0.5
        or dive_pullup_enable > 0.5
    )
    pullup_enable = float(point.family == "pullup" and point.load_factor > 1.0 + 1e-12)
    roll_enable = float(point.family == "roll" and abs(point.bank_angle_deg) > 1e-12)
    # For a constant-speed pull-up, n ~= 1 + V*q/g.  This is a command, not a
    # prescribed aircraft motion: MBDyn still solves the full rigid/flexible
    # dynamics and the achieved load factor is measured from the result.
    pullup_q_target = pullup_enable * (point.load_factor - 1.0) * 9.81 / point.velocity_mps
    # Feed-forward magnitude comes from the corrected SOL 144 trim map.  Its
    # sign is mapped to the MBDyn symmetric wing-flap convention here.
    bflap_1g_deg = float(np.interp(
        point.velocity_mps,
        [57.5, 65.0, 70.0],
        [0.8922637366, 0.7524614616, 0.6850896463],
    ))
    flap_increment = (
        -(point.load_factor - 1.0) * bflap_1g_deg * math.pi / 180.0
        if point.family == "pullup" else 0.0
    )
    text = text.replace(
        "set: const real V_INF =",
        (
            f"# MANEUVER_METADATA {json.dumps(asdict(point), sort_keys=True)}\n"
            f"# MANEUVER_MODEL_REVISION {model_revision(point)}\n"
            f"set: const real MANEUVER_LOAD_FACTOR = {point.load_factor:.10e};\n"
            f"set: const real MANEUVER_BANK_ANGLE = {point.bank_angle_deg:.10e}*deg2rad;\n"
            f"set: const real MANEUVER_ENABLE = {maneuver_enable:.1f};\n"
            f"set: const real PULLUP_ENABLE = {pullup_enable:.1f};\n"
            f"set: const real PULLUP_ANGLE_ENABLE = {angle_enable:.1f};\n"
            f"set: const real DIVE_PULLUP_ENABLE = {dive_pullup_enable:.1f};\n"
            f"set: const real ROLL_ENABLE = {roll_enable:.1f};\n"
            f"set: const real PULLUP_Q_TARGET = {pullup_q_target:.10e};\n"
            f"set: const real PULLUP_PITCH_ANGLE = {point.pitch_angle_deg:.10e}*deg2rad;\n"
            f"set: const real DIVE_PITCH_RATE_COMMAND = {point.pitch_rate_deg_s:.10e}*deg2rad;\n"
            f"set: const real DIVE_NOMINAL_LOAD_CLASS = {point.nominal_load_factor:.10e};\n"
            f"set: const real MANEUVER_FLAP_INCREMENT = {flap_increment:.10e};\n"
            f"set: const real MANEUVER_START = {maneuver_start:.10e};\n"
            f"set: const real MANEUVER_RISE_TIME = {maneuver_rise:.10e};\n"
            f"set: const real MANEUVER_RETURN_START = {maneuver_return:.10e};\n"
            f"set: const real MANEUVER_FALL_TIME = {maneuver_fall:.10e};\n"
            f"set: const real DIVE_END = {dive_end:.10e};\n"
            f"set: const real DIVE_PULLUP_START = {pullup_start:.10e};\n"
            f"set: const real DIVE_PULLUP_END = {pullup_end:.10e};\n"
            f"set: const real DIVE_RECOVERY_END = {recovery_end:.10e};\n"
            f"set: const real DIVE_MANEUVER_END = {dive_maneuver_end:.10e};\n"
            f"set: const real DIVE_CONTROL_RECOVERY_END = {dive_control_recovery_end:.10e};\n"
            "set: const integer MANEUVER_PROFILE_DRIVE = 12001;\n"
            "set: const integer MANEUVER_PITCH_COMMAND_DRIVE = 12002;\n"
            "set: const integer MANEUVER_Q_COMMAND_DRIVE = 12003;\n"
            "set: const integer MANEUVER_ROLL_COMMAND_DRIVE = 12004;\n"
            "set: const integer MANEUVER_FLAP_DRIVE = 12005;\n"
            "set: const integer MANEUVER_ALTITUDE_GATE_DRIVE = 12006;\n"
            "set: const integer MANEUVER_ANGLE_PROFILE_DRIVE = 12007;\n"
            "set: const integer MANEUVER_VZ_GATE_DRIVE = 12008;\n"
            "set: const integer DIVE_PULLUP_PROFILE_DRIVE = 12009;\n"
            "set: const integer DIVE_PULLUP_RATE_DRIVE = 12010;\n"
            "set: const integer DIVE_TRAJECTORY_GATE_DRIVE = 12011;\n"
            "set: const real V_INF ="
        ),
        1,
    )

    maneuver_drives = """
    # Actual time-domain maneuver commands.  Gravity remains physical (1 g).
    drive caller: MANEUVER_PROFILE_DRIVE,
        ramp, 0., 1./MANEUVER_RISE_TIME, MANEUVER_START, 1.;
    drive caller: MANEUVER_ANGLE_PROFILE_DRIVE,
        double ramp,
            1./MANEUVER_RISE_TIME,
            MANEUVER_START, MANEUVER_START + MANEUVER_RISE_TIME,
            -1./MANEUVER_FALL_TIME,
            MANEUVER_RETURN_START, MANEUVER_RETURN_START + MANEUVER_FALL_TIME,
            0.;
    # Normalized balanced trajectory: 0 -> -1 (dive) -> +1 (recovery) -> 0.
    drive caller: DIVE_PULLUP_PROFILE_DRIVE, string,
        "-((Time-MANEUVER_START)/(DIVE_END-MANEUVER_START))*((Time>=MANEUVER_START)&&(Time<DIVE_END))-((Time>=DIVE_END)&&(Time<DIVE_PULLUP_START))+(-1.+2.*(Time-DIVE_PULLUP_START)/(DIVE_PULLUP_END-DIVE_PULLUP_START))*((Time>=DIVE_PULLUP_START)&&(Time<DIVE_PULLUP_END))+((Time>=DIVE_PULLUP_END)&&(Time<DIVE_RECOVERY_END))+(1.-(Time-DIVE_RECOVERY_END)/(DIVE_MANEUVER_END-DIVE_RECOVERY_END))*((Time>=DIVE_RECOVERY_END)&&(Time<DIVE_MANEUVER_END))";
    drive caller: DIVE_PULLUP_RATE_DRIVE, string,
        "-((Time>=MANEUVER_START)&&(Time<DIVE_END))/(DIVE_END-MANEUVER_START)+2.*((Time>=DIVE_PULLUP_START)&&(Time<DIVE_PULLUP_END))/(DIVE_PULLUP_END-DIVE_PULLUP_START)-((Time>=DIVE_RECOVERY_END)&&(Time<DIVE_MANEUVER_END))/(DIVE_MANEUVER_END-DIVE_RECOVERY_END)";
    # Re-engage altitude/Vz smoothly after the attitude command has returned
    # to trim; a hard 0->1 switch creates a spurious rigid-body transient.
    drive caller: DIVE_TRAJECTORY_GATE_DRIVE, string,
        "((Time>=DIVE_MANEUVER_END)&&(Time<DIVE_CONTROL_RECOVERY_END))*(Time-DIVE_MANEUVER_END)/(DIVE_CONTROL_RECOVERY_END-DIVE_MANEUVER_END)+(Time>=DIVE_CONTROL_RECOVERY_END)";
    drive caller: MANEUVER_Q_COMMAND_DRIVE, string,
        "PULLUP_ENABLE*PULLUP_Q_TARGET*model::drive(MANEUVER_PROFILE_DRIVE,Time)+PULLUP_ANGLE_ENABLE*PULLUP_PITCH_ANGLE*(1./MANEUVER_RISE_TIME)*((Time>=MANEUVER_START)&&(Time<MANEUVER_START+MANEUVER_RISE_TIME))-PULLUP_ANGLE_ENABLE*PULLUP_PITCH_ANGLE*(1./MANEUVER_FALL_TIME)*((Time>=MANEUVER_RETURN_START)&&(Time<MANEUVER_RETURN_START+MANEUVER_FALL_TIME))+DIVE_PULLUP_ENABLE*PULLUP_PITCH_ANGLE*model::drive(DIVE_PULLUP_RATE_DRIVE,Time)";
    drive caller: MANEUVER_PITCH_COMMAND_DRIVE, string,
        "TRIM_PITCH+PULLUP_ENABLE*PULLUP_Q_TARGET*((Time>MANEUVER_START)*(Time-MANEUVER_START)-0.5*(Time>MANEUVER_START+MANEUVER_RISE_TIME)*MANEUVER_RISE_TIME-0.5*((Time>MANEUVER_START)&&(Time<MANEUVER_START+MANEUVER_RISE_TIME))*(Time-MANEUVER_START)*(2.-(Time-MANEUVER_START)/MANEUVER_RISE_TIME))+PULLUP_ANGLE_ENABLE*PULLUP_PITCH_ANGLE*model::drive(MANEUVER_ANGLE_PROFILE_DRIVE,Time)+DIVE_PULLUP_ENABLE*PULLUP_PITCH_ANGLE*model::drive(DIVE_PULLUP_PROFILE_DRIVE,Time)";
    drive caller: MANEUVER_ROLL_COMMAND_DRIVE, mult,
        const, ROLL_ENABLE*MANEUVER_BANK_ANGLE, reference, MANEUVER_PROFILE_DRIVE;
    drive caller: MANEUVER_FLAP_DRIVE, mult,
        const, MANEUVER_ENABLE*MANEUVER_FLAP_INCREMENT,
        reference, MANEUVER_PROFILE_DRIVE;
    drive caller: MANEUVER_ALTITUDE_GATE_DRIVE, string,
        "(1.-PULLUP_ENABLE*model::drive(MANEUVER_PROFILE_DRIVE,Time))*(1.-PULLUP_ANGLE_ENABLE*(Time>=MANEUVER_START))*(1.-DIVE_PULLUP_ENABLE+DIVE_PULLUP_ENABLE*model::drive(DIVE_TRAJECTORY_GATE_DRIVE,Time))";
    drive caller: MANEUVER_VZ_GATE_DRIVE, string,
        "(1.-max(PULLUP_ENABLE*model::drive(MANEUVER_PROFILE_DRIVE,Time),PULLUP_ANGLE_ENABLE*model::drive(MANEUVER_ANGLE_PROFILE_DRIVE,Time)))*(1.-DIVE_PULLUP_ENABLE+DIVE_PULLUP_ENABLE*model::drive(DIVE_TRAJECTORY_GATE_DRIVE,Time))";

"""
    marker = "    user defined: ALT_PID, pid,"
    if text.count(marker) != 1:
        raise RuntimeError("punto di inserimento dei drive di manovra non univoco")
    text = text.replace(marker, maneuver_drives + marker)
    text = _replace_pid_setpoint(text, "PITCH_PID", "reference, MANEUVER_PITCH_COMMAND_DRIVE,")
    text = _replace_pid_setpoint(text, "Q_PID", "reference, MANEUVER_Q_COMMAND_DRIVE,")
    text = _replace_pid_setpoint(text, "ROLL_PID", "reference, MANEUVER_ROLL_COMMAND_DRIVE,")

    lift_pattern = (
        r"drive caller:\s*LIFT_RAW_DRIVE, mult, const, LIFT_HOLD_GAIN,\s*"
        r"array, 2, reference, ALT_DRIVE, reference, VZ_DRIVE;"
    )
    lift_replacement = (
        "drive caller: LIFT_RAW_DRIVE, array, 3,\n"
        "        mult, const, LIFT_HOLD_GAIN,\n"
        "            mult, reference, MANEUVER_ALTITUDE_GATE_DRIVE, reference, ALT_DRIVE,\n"
        "        mult, const, LIFT_HOLD_GAIN,\n"
        "            mult, reference, MANEUVER_VZ_GATE_DRIVE, reference, VZ_DRIVE,\n"
        "        reference, MANEUVER_FLAP_DRIVE;"
    )
    text, count = re.subn(lift_pattern, lift_replacement, text)
    if count != 1:
        raise RuntimeError("LIFT_RAW_DRIVE baseline non trovato")
    text = _replace_real(text, "BFF_RAP_AMPLITUDE", 0.20 * math.pi / 180.0 if point.excited else 0.0)

    # The validated controller remains active before release.  During SAS-off
    # every aerodynamic surface, including both body flaps, is frozen at its
    # value at release.  The aircraft therefore evolves open-loop from the
    # same maneuver state in shadow and excited runs, with no pilot/SAS command
    # capable of contaminating the BFF band.
    surface_include = f'include: "{BASELINE / "INCLUDE/control_surfaces.mbd"}";'
    surfaces = (BASELINE / "INCLUDE/control_surfaces.mbd").read_text()
    if text.count(surface_include) != 1:
        raise RuntimeError("include delle superfici non univoco")
    text = text.replace(
        surface_include,
        "# MANEUVER_CONTROL: all aerodynamic surfaces held at release\n" + surfaces,
    )

    required = (
        "gravity: uniform, 0., 0., -1., const, GRAVITY",
        "reference, MANEUVER_PITCH_COMMAND_DRIVE",
        "reference, MANEUVER_Q_COMMAND_DRIVE",
        "reference, MANEUVER_ROLL_COMMAND_DRIVE",
        "reference, MANEUVER_ALTITUDE_GATE_DRIVE",
        "reference, MANEUVER_VZ_GATE_DRIVE",
        "MANEUVER_CONTROL: all aerodynamic surfaces held at release",
        "reference, BFL_HOLD_DRIVE",
        "reference, BFR_HOLD_DRIVE",
        "node, BASE_NODE, structural, string, \"Omega[2]\", direct",
        "force: NASTRAN_DLM_ROM_FORCE, modal, MODAL_JOINT",
        "force: SAS_MODAL_DAMPER_FORCE, modal, MODAL_JOINT",
        "SAS_OFF_DURATION = 2.0500000000e+00",
    )
    if not all(item in text for item in required):
        raise RuntimeError("invarianti del caso maneuver-BFF non soddisfatte")
    return text


def run_point(point: ManeuverPoint, output: Path, overwrite: bool = False) -> Path:
    # MBDyn is launched with cwd=BASELINE so all user-provided relative output
    # paths must be resolved before building the input and output arguments.
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    stem = point.stem
    prefix = output / stem
    mbd = prefix.with_suffix(".mbd")
    nc = prefix.with_suffix(".nc")
    rendered = render_point(point)
    if nc.exists() and not overwrite:
        raise FileExistsError(f"{nc} esiste; usare --overwrite")
    if overwrite:
        for suffix in (".mbd", ".nc", ".out", ".log", ".stdout", ".mov"):
            path = prefix.with_suffix(suffix)
            if path.is_file():
                path.unlink()
    mbd.write_text(rendered)
    result = subprocess.run(
        [str(MBDYN), "-s", "-f", str(mbd), "-o", str(prefix)],
        cwd=BASELINE,
        text=True,
        capture_output=True,
    )
    prefix.with_suffix(".stdout").write_text(result.stdout + result.stderr)
    if result.returncode or not nc.exists():
        raise RuntimeError(f"MBDyn fallito per {stem}; vedere {prefix.with_suffix('.stdout')}")
    if (
        point.load_factor == 1.0
        and abs(point.bank_angle_deg) < 1e-12
        and abs(point.pitch_angle_deg) < 1e-12
    ):
        baseline_module().validate_release_state(nc, rendered)
    else:
        # A maneuver must not satisfy the level-flight q/Vz limits.  Validate
        # the actual body-normal load reached at release instead.
        from analyse_time_domain_pairs import trajectory_metrics
        audit = trajectory_metrics(nc, mbd)
        load_error = audit["achieved_n_mean_sas_off"] - point.load_factor
        if point.family in {"pullup_angle", "dive_pullup"}:
            print(
                f"[maneuver {point.family}] finestra SAS-off: "
                f"ampiezza={point.pitch_angle_deg:.3f} deg, "
                f"q_command={point.pitch_rate_deg_s:.3f} deg/s, "
                f"n={audit['achieved_n_mean_sas_off']:.4f}, "
                f"q={audit['q_release_deg_s']:+.4f} deg/s"
            )
        elif abs(load_error) > 0.075:
            print(
                "[maneuver WARNING] punto escluso dall'identificazione: "
                f"n medio SAS-off={audit['achieved_n_mean_sas_off']:.4f}, "
                f"target={point.load_factor:.4f}",
                file=sys.stderr,
            )
        else:
            print(
                "[maneuver] finestra SAS-off valida: "
                f"n={audit['achieved_n_mean_sas_off']:.4f}, "
                f"q={audit['q_release_deg_s']:+.4f} deg/s, "
                f"p={audit['p_release_deg_s']:+.4f} deg/s"
            )
    return nc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=("pullup", "pullup_angle", "dive_pullup", "roll"), required=True)
    parser.add_argument("--velocity", type=float, required=True)
    parser.add_argument("--load-factor", type=float, required=True)
    parser.add_argument("--bank-angle", type=float, default=0.0)
    parser.add_argument("--pitch-angle", type=float, default=0.0,
                        help="incremento/ampiezza di pitch [deg] per pullup_angle o dive_pullup")
    parser.add_argument("--pitch-rate", type=float, default=0.0,
                        help="pitch-rate nominale [deg/s] per dive_pullup")
    parser.add_argument("--nominal-load", type=float, default=1.0,
                        help="sola etichetta della classe di aggressivita'; n e' misurato a posteriori")
    parser.add_argument("--shadow", action="store_true", help="stessa manovra senza rap BFF")
    parser.add_argument("--output", type=Path, default=Path(CONFIG["default_output"]) / "single_cases")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    point = ManeuverPoint(
        args.family, args.velocity, args.load_factor, args.bank_angle,
        not args.shadow, args.pitch_angle, args.pitch_rate, args.nominal_load,
    )
    if args.dry_run:
        output = args.output.expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        path = output / f"{point.stem}.mbd"
        path.write_text(render_point(point))
        print(path)
        return
    print(run_point(point, args.output, args.overwrite))


if __name__ == "__main__":
    main()
