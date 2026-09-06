#!/usr/bin/env python3
"""Case construction shared by the MANOUVER_STIFNESS tools."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONFIG = json.loads((ROOT / "campaign_config.json").read_text())
MODEL = CONFIG["model"]
def _model_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (ROOT / path).resolve()


MANEUVER_DIR = _model_path(MODEL["maneuver_directory"])
BASELINE_DIR = _model_path(MODEL["baseline_directory"])
MODEL_REVISION = "manouver_stifness_v2_reduced_window_and_output"


def _load_module(name: str, path: Path):
    if str(path.parent) not in sys.path:
        sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"impossibile caricare {path}")
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolves postponed annotations through sys.modules while
    # the module body is executing.
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


MANEUVER = _load_module("stiffness_source_maneuver_case", MANEUVER_DIR / "maneuver_case.py")


@dataclass(frozen=True)
class StudyCase:
    campaign: str
    velocity_mps: float
    nominal_load_factor: float
    time_step_s: float
    mode7_frequency_scale: float
    excited: bool

    @property
    def maneuver_point(self):
        q_deg_s, amplitude_deg = command_for_class(
            self.velocity_mps, self.nominal_load_factor
        )
        return MANEUVER.ManeuverPoint(
            "dive_pullup",
            self.velocity_mps,
            1.0,
            0.0,
            self.excited,
            amplitude_deg,
            q_deg_s,
            self.nominal_load_factor,
        )

    @property
    def stem(self) -> str:
        dt = tag(self.time_step_s, 4)
        scale = tag(self.mode7_frequency_scale, 3)
        return f"{self.campaign}__dt_{dt}__f7_{scale}__{self.maneuver_point.stem}"

    def metadata(self) -> dict:
        data = asdict(self)
        q_deg_s, amplitude_deg = command_for_class(
            self.velocity_mps, self.nominal_load_factor
        )
        data.update(
            model_revision=MODEL_REVISION,
            pitch_rate_command_deg_s=q_deg_s,
            pitch_amplitude_deg=amplitude_deg,
            stiffness_delta_k7=stiffness_delta_k7(self.mode7_frequency_scale),
            modification=(
                "none_linear_modal_stiffness"
                if abs(self.mode7_frequency_scale - 1.0) < 1e-12
                else "parametric_mode7_stiffness_force_not_physical_prestress"
            ),
        )
        return data


def tag(value: float, digits: int) -> str:
    return f"{value:.{digits}f}".replace("-", "m").replace(".", "p")


def command_for_class(velocity_mps: float, nominal_load: float) -> tuple[float, float]:
    """Return pitch-rate command and centered dive amplitude.

    n is only a design label.  The achieved value is always reconstructed
    from the shadow trajectory.
    """
    q_deg_s = (
        (nominal_load - 1.0)
        * float(CONFIG["physics"]["gravity_m_s2"])
        / velocity_mps
        * 180.0
        / math.pi
    )
    source_timing = MANEUVER.CONFIG["dive_pullup_timing"]
    centered_duration = float(source_timing["pullup_release_lead_s"]) + float(
        source_timing["sas_off_duration_s"]
    )
    return q_deg_s, 0.5 * q_deg_s * centered_duration


def stiffness_delta_k7(frequency_scale: float) -> float:
    """Added modal stiffness for a requested natural-frequency scale.

    Modal mass is one in the current FEM.  The injected generalized force is
    -DeltaK*q7, hence DeltaK=(s_f**2-1)*K77; K77 is read from RECORD GROUP 10.
    """
    k77 = float(MODEL["dominant_mode_modal_stiffness"])
    return (frequency_scale * frequency_scale - 1.0) * k77


def first_order_lowpass(cutoff_hz: float, dt: float) -> tuple[float, float, float]:
    # Prewarp the declared digital cutoff before the bilinear transform.  This
    # reproduces the validated dt=0.01 coefficients exactly and preserves the
    # same -3 dB frequency when dt changes.
    omega = (2.0 / dt) * math.tan(math.pi * cutoff_hz * dt)
    denominator = 2.0 + omega * dt
    return (
        (2.0 - omega * dt) / denominator,
        omega * dt / denominator,
        omega * dt / denominator,
    )


def second_order_butterworth(cutoff_hz: float, dt: float) -> tuple[float, ...]:
    """Bilinear discretization in MBDyn's discrete-filter sign convention."""
    omega = (2.0 / dt) * math.tan(math.pi * cutoff_hz * dt)
    rate = 2.0 / dt
    a0 = rate * rate + math.sqrt(2.0) * omega * rate + omega * omega
    a1 = -2.0 * rate * rate + 2.0 * omega * omega
    a2 = rate * rate - math.sqrt(2.0) * omega * rate + omega * omega
    b0 = omega * omega / a0
    return (-a1 / a0, -a2 / a0, b0, 2.0 * b0, b0)


def actuator_coefficients(tau: float, dt: float) -> tuple[float, float, float]:
    denominator = 2.0 * tau + dt
    return (
        (2.0 * tau - dt) / denominator,
        dt / denominator,
        dt / denominator,
    )


def _replace_real(text: str, name: str, value: float) -> str:
    pattern = rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*[^;]+;"
    text, count = re.subn(pattern, f"set: const real {name} = {value:.15e};", text)
    if count != 1:
        raise RuntimeError(f"{name}: attesa una definizione, trovate {count}")
    return text


def _rediscretize(text: str, dt: float) -> str:
    digital = CONFIG["digital_systems"]
    lp = second_order_butterworth(float(digital["second_order_lowpass_hz"]), dt)
    alt = first_order_lowpass(float(digital["altitude_lowpass_hz"]), dt)
    vz = first_order_lowpass(float(digital["vertical_speed_lowpass_hz"]), dt)
    pitch_rate = first_order_lowpass(float(digital["pitch_rate_lowpass_hz"]), dt)
    actuator = actuator_coefficients(float(digital["actuator_time_constant_s"]), dt)
    coefficients = {
        "LP_A1": lp[0], "LP_A2": lp[1], "LP_B0": lp[2],
        "LP_B1": lp[3], "LP_B2": lp[4],
        "ALT_LP_A1": alt[0], "ALT_LP_B0": alt[1], "ALT_LP_B1": alt[2],
        "VZ_LP_A1": vz[0], "VZ_LP_B0": vz[1], "VZ_LP_B1": vz[2],
        "Q_LP_A1": pitch_rate[0], "Q_LP_B0": pitch_rate[1],
        "Q_LP_B1": pitch_rate[2],
        "ACT_A1": actuator[0], "ACT_B0": actuator[1], "ACT_B1": actuator[2],
    }
    text = _replace_real(text, "TIME_STEP", dt)
    for name, value in coefficients.items():
        text = _replace_real(text, name, value)
    return text


def _inject_stiffness_screen(text: str, case: StudyCase) -> str:
    if abs(case.mode7_frequency_scale - 1.0) < 1e-12:
        return text
    constants_marker = "set: const integer MODAL_JOINT = 5;"
    if text.count(constants_marker) != 1:
        raise RuntimeError("costante MODAL_JOINT non trovata in modo univoco")
    constants = f"""

# Parametric stiffness screen.  This is deliberately not called prestress:
# it changes only K77 and is used as a decision gate for a future physical ROM.
set: const integer STIFFNESS_SCREEN_FORCE = 9601;
set: const real STIFFNESS_SCREEN_FREQUENCY_SCALE = {case.mode7_frequency_scale:.15e};
set: const real STIFFNESS_SCREEN_DELTA_K7 = {stiffness_delta_k7(case.mode7_frequency_scale):.15e};
set: const real STIFFNESS_SCREEN_RAMP_TIME = 5.;
"""
    text = text.replace(constants_marker, constants_marker + constants)
    text, count = re.subn(
        r"(?m)^(\s*forces:\s*)2(\s*;)", r"\g<1>3\g<2>", text, count=1
    )
    if count != 1:
        raise RuntimeError("conteggio delle forze baseline non trovato")
    marker = "    # Idealized active modal damper"
    if text.count(marker) != 1:
        raise RuntimeError("punto di inserimento della forza modale non univoco")
    force = """    # Generalized force -DeltaK*q7.  The five-second ramp lets the
    # controlled aircraft settle to the modified stiffness before the maneuver.
    force: STIFFNESS_SCREEN_FORCE, modal, MODAL_JOINT,
        list, 1, 7,
        array, 1,
            element, MODAL_JOINT, joint, string, "q[7]",
                string, "-min(1.,Time/STIFFNESS_SCREEN_RAMP_TIME)*STIFFNESS_SCREEN_DELTA_K7*Var",
        output, no;

"""
    return text.replace(marker, force + marker)


def _constant(text: str, name: str) -> float:
    match = re.search(
        rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*([-+0-9.eE]+)\s*;",
        text,
    )
    if not match:
        raise RuntimeError(f"{name}: costante numerica non trovata")
    return float(match.group(1))


def _optimize_production_output(text: str, case: StudyCase) -> str:
    """Remove integration and output that cannot affect the identified pole."""
    execution = CONFIG["execution"]
    if execution["final_time_policy"] != "identification_window_only":
        raise RuntimeError("final_time_policy non supportata")
    release = _constant(text, "SAS_OFF_START")
    duration = _constant(text, "SAS_OFF_DURATION")
    margin = int(execution["post_sas_off_margin_steps"]) * case.time_step_s
    text = _replace_real(text, "FINAL_TIME", release + duration + margin)

    constants_marker = "set: const integer MODAL_JOINT = 5;"
    output_constants = (
        f"\nset: const real STUDY_OUTPUT_PREHISTORY = "
        f"{float(execution['output_prehistory_s']):.15e};\n"
        f"set: const real STUDY_OUTPUT_SAMPLE_STEP = "
        f"{float(execution['output_sample_step_s']):.15e};\n"
    )
    if text.count(constants_marker) != 1:
        raise RuntimeError("costante MODAL_JOINT non univoca")
    text = text.replace(constants_marker, constants_marker + output_constants, 1)

    control_marker = "    print: dof stats;"
    output_control = (
        "    # Production output: fixed physical sampling, only near release.\n"
        "    default output: none;\n"
        "    output meter: closest next, SAS_OFF_START - STUDY_OUTPUT_PREHISTORY, "
        "forever, const, STUDY_OUTPUT_SAMPLE_STEP;\n"
    )
    if text.count(control_marker) != 1:
        raise RuntimeError("control data marker non univoco")
    text = text.replace(control_marker, output_control + control_marker, 1)

    node_include = f'    include: "{BASELINE_DIR / "INCLUDE/node.mbd"}";'
    nodes_marker = node_include + "\nend: nodes;"
    node_outputs = "".join(
        f"    output: structural, {int(label)};\n"
        for label in execution["output_structural_nodes"]
    )
    if text.count(nodes_marker) != 1:
        raise RuntimeError("sezione nodi non trovata in modo univoco")
    text = text.replace(
        nodes_marker,
        node_include + "\n" + node_outputs + "end: nodes;",
        1,
    )

    aero_include = f'    include: "{BASELINE_DIR / "INCLUDE/aerobody.mbd"}";'
    elements_marker = aero_include + "\nend: elements;"
    joint_outputs = "".join(
        f"    output: joint, {int(label)};\n"
        for label in execution["output_joints"]
    )
    if text.count(elements_marker) != 1:
        raise RuntimeError("fine sezione elementi non trovata in modo univoco")
    text = text.replace(
        elements_marker,
        aero_include + "\n" + joint_outputs + "end: elements;",
        1,
    )
    return text


def render_case(case: StudyCase) -> str:
    if case.campaign not in CONFIG["campaigns"]:
        raise ValueError(f"campagna sconosciuta: {case.campaign}")
    point = case.maneuver_point
    text = MANEUVER.render_point(point)
    text = _rediscretize(text, case.time_step_s)
    text = _inject_stiffness_screen(text, case)
    text = _optimize_production_output(text, case)
    metadata = json.dumps(case.metadata(), sort_keys=True, separators=(",", ":"))
    text = text.replace(
        "# MANEUVER_METADATA ",
        f"# STIFFNESS_STUDY_METADATA {metadata}\n# MANEUVER_METADATA ",
        1,
    )
    required = (
        f"TIME_STEP = {case.time_step_s:.15e}",
        "MANEUVER_CONTROL: all aerodynamic surfaces held at release",
        "force: NASTRAN_DLM_ROM_FORCE, modal, MODAL_JOINT",
        "force: SAS_MODAL_DAMPER_FORCE, modal, MODAL_JOINT",
        "STIFFNESS_STUDY_METADATA",
        "default output: none",
        "output meter: closest next",
        "output: joint, 5",
    )
    if not all(item in text for item in required):
        raise RuntimeError(f"invarianti del caso non soddisfatte per {case.stem}")
    return text


def build_cases(campaign_name: str) -> list[StudyCase]:
    cfg = CONFIG["campaigns"].get(campaign_name)
    if cfg is None:
        raise ValueError(f"campagna sconosciuta: {campaign_name}")
    cases = [
        StudyCase(campaign_name, float(v), float(n), float(dt), float(scale), excited)
        for v in cfg["velocities_mps"]
        for n in cfg["nominal_load_factors"]
        for dt in cfg["time_steps_s"]
        for scale in cfg["mode7_frequency_scales"]
        for excited in (False, True)
    ]
    stems = [case.stem for case in cases]
    if len(stems) != len(set(stems)):
        raise RuntimeError("stem duplicati nella matrice dei casi")
    return cases


def source_fingerprints() -> dict[str, str]:
    paths = {
        "campaign_generator": Path(__file__).resolve(),
        "campaign_config": ROOT / "campaign_config.json",
        "maneuver_case": MANEUVER_DIR / "maneuver_case.py",
        "study_config": MANEUVER_DIR / "study_config.json",
        "baseline_run_case": BASELINE_DIR / "run_case.py",
        "baseline_main": BASELINE_DIR / "main_bff_open_loop.mbd",
        "fem": Path(MODEL["fem_file"]),
    }
    return {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in paths.items()
    }
