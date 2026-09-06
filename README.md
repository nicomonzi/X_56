# X-56 aeroelastic research workspace

Repository riorganizzata a partire da `TESI` il 5 settembre 2026. La sorgente
originale non è stata modificata. Il rapporto storico completo è in
`docs/report/RAPPORTO_GENERALE.md`; le decisioni di migrazione sono in
`docs/decision_log/MIGRAZIONE_DA_TESI.md`.

## Workflow attivi

- `workflows/bff_open_loop`: baseline BFF steady con SAS-off;
- `workflows/maneuver_bff`: envelope, campagne paired e studio prestress;
- `workflows/trim`: trim longitudinale MBDyn;
- `workflows/coupled_dust`: coupling MBDyn–DUST e staging trim 55 m/s.

I modelli e i dati condivisi sono sotto `models/` e `assets/`; le verifiche
scientifiche sono sotto `validation/`. `archive/` non è una sorgente attiva.

## Verifica dopo il clone

```bash
python3 tools/verify_dependencies.py
python3 workflows/bff_open_loop/run_case.py --help
python3 workflows/maneuver_bff/run_sweep.py --help
python3 workflows/trim/run_trim_sweep.py --help
python3 workflows/coupled_dust/run_case.py --help
```

Le simulazioni non partono con questi comandi. DUST/preCICE richiedono inoltre
binari e binding configurati in `workflows/coupled_dust/config/machine.env`.
Gli output pesanti devono restare fuori Git e vanno associati a manifest e hash.
