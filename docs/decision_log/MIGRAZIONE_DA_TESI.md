# Migrazione da TESI

Data: 5 settembre 2026. Sorgente: `/home/nicomonzi/TESI`. Destinazione: `/home/nicomonzi/X_56`.

## Decisioni

- Il contenuto della sorgente è stato copiato, mai spostato.
- `BFF_maneuver_envelope` e `MANOUVER_STIFNESS` sono stati uniti in
  `workflows/maneuver_bff` perché il secondo importa direttamente il primo.
- I FEM identici sono conservati una sola volta in `assets/fem`; i vecchi nomi
  sono link simbolici relativi, compatibili con Linux/WSL.
- I BULK R11 canonici sono in `assets/nastran_bulk/x56_r11`; i casi Nastran
  attivi mantengono il nome `BULK` tramite link relativo.
- Cache, virtualenv, build, tentativi incompleti e output grezzi NC/AER/MOV/MOD/
  JNT/USR non sono stati copiati. Manifest, CSV, JSON, report, figure e input
  sorgente sono stati conservati.
- `archive` è genealogia scientifica: non è incluso nel gate di dipendenze dei
  workflow attivi.
- I percorsi esterni ZENO e Desktop restano configurabili perché contengono
  risultati non inclusi nel repository; non sono dipendenze sorgente interne.

## Mappa principale

| TESI | X_56 |
|---|---|
| `BFF_open_loop` | `workflows/bff_open_loop` |
| `BFF_maneuver_envelope` + `MANOUVER_STIFNESS` | `workflows/maneuver_bff` |
| `TRIM` | `workflows/trim` |
| `BFF_DUST_55` | `workflows/coupled_dust` |
| `DUST` | `models/dust/x56` |
| `NASTRAN/REALASED_MODEL` | `models/nastran/x56_r11` |
| `NASTRAN_SIMULATIONS` | `validation/modal` e `validation/static_5g` |
| `X56_AERO_POLAR` | `validation/aero_polar` |
| `MBDYN`, `TEST`, vecchi `bff_*` | `archive` |

## Limiti mantenuti espliciti

La produzione DUST resta bloccata; la validazione statica 5 g non soddisfa la
soglia esterna; il gate fisico del prestress è rettificato e i controlli causali
sono ancora necessari. La riorganizzazione non cambia questi verdetti.
