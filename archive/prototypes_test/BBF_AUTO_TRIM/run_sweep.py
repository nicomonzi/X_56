#!/usr/bin/env python3
"""Esegue in sequenza i casi MBDyn da 23 a 36 m/s, senza analisi."""

from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent
MAIN = ROOT / "main_bbf.mbd"
MBDYN = "/usr/local/mbdyn/bin/mbdyn"
OUTPUT = Path("/mnt/c/Users/Utente/Desktop/RESULTS")


def main() -> None:
    source = MAIN.read_text(encoding="utf-8")
    pattern = r"(?m)^set:\s*const\s+real\s+V_INF\s*=\s*[^;]+;"

    for velocity in range(23, 37):
        case = OUTPUT / f"V_{velocity:03d}_mps"
        case.mkdir(parents=True, exist_ok=True)

        # Ogni caso conserva l'input esatto usato, senza modificare main_bbf.mbd.
        case_source, replacements = re.subn(
            pattern,
            f"set: const real V_INF = {velocity:.1f};",
            source,
            count=1,
        )
        if replacements != 1:
            raise RuntimeError("Impossibile trovare un'unica definizione di V_INF")

        case_input = case / "case_input.mbd"
        case_input.write_text(case_source, encoding="utf-8")
        output_prefix = case / "case"

        # MBDyn risolve "./INCLUDE" rispetto alla cartella del file di input.
        # La copia temporanea resta quindi accanto a main_bbf.mbd, mentre tutti
        # i risultati e la copia permanente dell'input vengono scritti su C:.
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".run_sweep_",
            suffix=".mbd",
            dir=ROOT,
            delete=False,
        ) as temporary:
            temporary.write(case_source)
            active_input = Path(temporary.name)

        print(f"\n=== V_INF = {velocity} m/s ===", flush=True)
        try:
            subprocess.run(
                [MBDYN, "-f", str(active_input), "-o", str(output_prefix)],
                cwd=ROOT,
                check=True,
            )
        finally:
            active_input.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
