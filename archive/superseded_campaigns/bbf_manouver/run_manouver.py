#!/usr/bin/env python3
"""Esegue un caso X-56: livellato, pull-up o roll, con burst BFF opzionale."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / "main_bbf.mbd"
OUTPUT = ROOT / "output"
MBDYN = Path("/usr/local/mbdyn/bin/mbdyn")


def replace_real(text: str, name: str, value: float) -> str:
    text, count = re.subn(
        rf"(?m)^set:\s*const\s+real\s+{re.escape(name)}\s*=\s*[^;]+;",
        f"set: const real {name} = {value:.8f};",
        text,
    )
    if count != 1:
        raise RuntimeError(f"{name}: attesa una definizione, trovate {count}")
    return text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("level", "pullup", "roll"), default="level")
    parser.add_argument("--velocity", type=float, default=60.8421, help="V_INF [m/s]")
    parser.add_argument(
        "--bff", action=argparse.BooleanOptionalAction, default=True,
        help="abilita/disabilita il burst simmetrico BFF",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="genera il caso senza eseguire MBDyn")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.velocity <= 0.0:
        raise ValueError("--velocity deve essere positiva")

    flags = {
        "PULLUP_ENABLE": float(args.mode == "pullup"),
        "ROLL_ENABLE": float(args.mode == "roll"),
        "BFF_ENABLE": float(args.bff),
    }
    text = replace_real(MODEL.read_text(), "V_INF", args.velocity)
    for name, value in flags.items():
        text = replace_real(text, name, value)
    # Il caso renderizzato vive in output/: gli include devono continuare a
    # puntare alla cartella del modello e non alla directory del file generato.
    text = text.replace('"./INCLUDE/', f'"{ROOT / "INCLUDE"}/')

    OUTPUT.mkdir(exist_ok=True)
    velocity_tag = f"{args.velocity:08.4f}".replace(".", "p")
    bff_tag = "bff" if args.bff else "no_bff"
    stem = f"{args.mode}_{bff_tag}_V_{velocity_tag}"
    input_path = OUTPUT / f"{stem}.mbd"
    prefix = OUTPUT / stem
    nc_path = OUTPUT / f"{stem}.nc"

    if nc_path.exists() and not args.overwrite:
        raise FileExistsError(f"{nc_path} esiste; usa --overwrite")
    if args.overwrite:
        for suffix in (".mbd", ".nc", ".out", ".log", ".stdout"):
            path = OUTPUT / f"{stem}{suffix}"
            if path.is_file():
                path.unlink()

    input_path.write_text(text)
    print(f"[case] mode={args.mode}, BFF={args.bff}, V={args.velocity:.4f} m/s")
    print(f"[input] {input_path}")
    if args.dry_run:
        return

    result = subprocess.run(
        [str(MBDYN), "-s", "-f", str(input_path), "-o", str(prefix)],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    (OUTPUT / f"{stem}.stdout").write_text(result.stdout + result.stderr)
    print(f"[done] returncode={result.returncode}, output={prefix}")
    if result.returncode:
        raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
