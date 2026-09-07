#!/usr/bin/env python3
"""Inject the generated prestress ROM into an already-rendered maneuver MBDyn case."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load_config(path: Path | None) -> dict:
    return json.loads((path or ROOT / "config.json").expanduser().resolve().read_text())


def nz_drive_text(config: dict, source: str) -> str:
    model = config["model"]
    drive = config["mbdyn"]["delta_n_drive"]
    if source == "commanded":
        return (
            "    # Scheduled Nastran load class; safest source for the first A/B campaign.\n"
            f"    drive caller: {drive}, string,\n"
            '        "(MANEUVER_LOAD_FACTOR-1.)*model::drive(MANEUVER_PROFILE_DRIVE,Time)";\n'
        )
    if source == "measured":
        return (
            "    # Instantaneous body-z specific load factor. This closes an acceleration-to-force loop.\n"
            "    drive caller: PRESTRESS_ROLL_DRIVE,\n"
            f'        node, {model["base_node"]}, structural, string, "E[1]", direct;\n'
            "    drive caller: PRESTRESS_PITCH_DRIVE,\n"
            f'        node, {model["base_node"]}, structural, string, "E[2]", direct;\n'
            "    drive caller: PRESTRESS_BODY_AZ_DRIVE,\n"
            f'        node, {model["base_node"]}, structural, string, "xPP[3]", direct;\n'
            f"    drive caller: {drive}, string,\n"
            '        "model::drive(PRESTRESS_BODY_AZ_DRIVE,Time)/GRAVITY+cos(model::drive(PRESTRESS_ROLL_DRIVE,Time))*cos(model::drive(PRESTRESS_PITCH_DRIVE,Time))-1.";\n'
        )
    if source == "maneuver_pitch_rate":
        return (
            "    # Exogenous load schedule used by the validated dive-pull-up campaign.\n"
            "    # q_cmd is in rad/s; VINF/GRAVITY converts it to Delta n.  The\n"
            "    # schedule follows dive, pull-up and SAS-off without closing an\n"
            "    # acceleration-to-generalized-force algebraic loop.\n"
            f"    drive caller: {drive}, string,\n"
            '        "DIVE_PULLUP_ENABLE*VINF*model::drive(MANEUVER_Q_COMMAND_DRIVE,Time)/GRAVITY";\n'
        )
    raise ValueError(
        "nz-source deve essere commanded, measured o maneuver_pitch_rate"
    )


def render_text(
    text: str,
    config_path: Path | None,
    source: str | None,
    generated_dir: Path | None = None,
) -> str:
    config = load_config(config_path)
    source = source or config["mbdyn"]["load_factor_source"]
    generated = (generated_dir or ROOT / "generated").expanduser().resolve()
    constants_path = generated / "prestress_constants.mbd"
    force_path = generated / "prestress_modal_force.mbd"
    if not constants_path.is_file() or not force_path.is_file():
        raise FileNotFoundError("prima eseguire: python3 prestress_rom.py")
    if "PRESTRESS_ROM_INJECTED" in text:
        raise ValueError("il caso contiene gia la ROM prestress")
    metadata_match = re.search(r"(?m)^# MANEUVER_METADATA (\{.*\})$", text)
    if source in {"commanded", "maneuver_pitch_rate"}:
        if metadata_match is None:
            raise ValueError(f"nz-source {source} richiede un caso maneuver_bff renderizzato")
        family = json.loads(metadata_match.group(1)).get("family")
        expected_family = "pullup" if source == "commanded" else "dive_pullup"
        if family != expected_family:
            raise ValueError(
                f"nz-source {source} richiede family={expected_family}, non {family!r}"
            )
    constants = constants_path.read_text()
    force = force_path.read_text()

    initial_marker = "begin: initial value;"
    if text.count(initial_marker) != 1:
        raise ValueError("marker initial value non univoco")
    text = text.replace(
        initial_marker,
        "# PRESTRESS_ROM_INJECTED\n"
        + f"# PRESTRESS_ROM_METADATA {json.dumps({'enabled': True, 'nz_source': source}, sort_keys=True)}\n"
        + constants + "\n" + initial_marker,
        1,
    )
    text, count = re.subn(
        r"(?m)^(\s*forces:\s*)(\d+)(\s*;)",
        lambda match: f"{match.group(1)}{int(match.group(2)) + 1}{match.group(3)}",
        text,
        count=1,
    )
    if count != 1:
        raise ValueError("conteggio forces non trovato")
    marker = "    # Idealized active modal damper"
    if text.count(marker) != 1:
        raise ValueError("marker della forza modale non univoco")
    insertion = nz_drive_text(config, source) + "\n" + "    " + force.replace("\n", "\n    ").rstrip() + "\n\n"
    text = text.replace(marker, insertion + marker, 1)
    required = (
        "PRESTRESS_ROM_INJECTED",
        f"force: {config['mbdyn']['force_label']}, modal",
        f"drive caller: {config['mbdyn']['delta_n_drive']}",
        "force: NASTRAN_DLM_ROM_FORCE, modal",
    )
    if not all(item in text for item in required):
        raise RuntimeError("invarianti di rendering non soddisfatte")
    return text


def render(
    input_path: Path,
    output_path: Path,
    config_path: Path | None,
    source: str | None,
    generated_dir: Path | None = None,
) -> None:
    text = render_text(
        input_path.expanduser().resolve().read_text(),
        config_path,
        source,
        generated_dir,
    )
    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="caso maneuver_bff gia renderizzato")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--generated", type=Path, help="directory prodotta da prestress_rom.py")
    parser.add_argument(
        "--nz-source", choices=("commanded", "measured", "maneuver_pitch_rate")
    )
    args = parser.parse_args()
    render(args.input, args.output, args.config, args.nz_source, args.generated)
    print(args.output.expanduser().resolve())


if __name__ == "__main__":
    main()
