# Rapporti del progetto TESI

Questo dossier fotografa la cartella `/home/nicomonzi/TESI` al 5 settembre
2026. Il punto di ingresso è [RAPPORTO_GENERALE.md](RAPPORTO_GENERALE.md).

Contenuto:

- `RAPPORTO_GENERALE.md`: stato del progetto, cronologia, cartelle obsolete,
  duplicati, rischi e riorganizzazione proposta;
- `METODOLOGIA.md`: criteri, definizioni e limiti dell'audit;
- `cartelle/`: un rapporto interpretativo per ogni sottocartella di primo
  livello del progetto;
- `dati/inventario_file.csv`: tutti i file regolari con dimensione, data,
  SHA-256 e stato di tracciamento Git;
- `dati/duplicati_sha256.csv`: catalogo completo delle copie byte-per-byte;
- `dati/riepilogo_cartelle.csv`: statistiche quantitative per cartella;
- `dati/inventario_sottocartelle.csv`: conteggi, dimensioni e tipi di file per
  ogni directory annidata, comprese quelle vuote;
- `dati/catalogo_codice_python.csv`: righe, docstring, classi, funzioni e import
  di ogni script Python del progetto, esclusi ambienti e build di terze parti;
- `dati/cronologia_git.csv`: i 68 commit dal più recente al più vecchio;
- `dati/snapshot.json`: stato macchina usato per confrontare gli aggiornamenti;
- `PROGRESSI.md`: registro append-only delle variazioni rilevate.

Rapporti delle cartelle:

- [`BFF_DUST_55`](cartelle/BFF_DUST_55.md)
- [`BFF_maneuver_envelope`](cartelle/BFF_maneuver_envelope.md)
- [`BFF_open_loop`](cartelle/BFF_open_loop.md)
- [`DUST`](cartelle/DUST.md)
- [`GRAFICI`](cartelle/GRAFICI.md)
- [`MANOUVER_STIFNESS`](cartelle/MANOUVER_STIFNESS.md)
- [`MBDYN`](cartelle/MBDYN.md)
- [`NASTRAN`](cartelle/NASTRAN.md)
- [`NASTRAN_SIMULATIONS`](cartelle/NASTRAN_SIMULATIONS.md)
- [`TEST`](cartelle/TEST.md)
- [`TRIM`](cartelle/TRIM.md)
- [`X56_AERO_POLAR`](cartelle/X56_AERO_POLAR.md)
- [`_trim_coupled_stage`](cartelle/_trim_coupled_stage.md)
- [`bbf_manouver`](cartelle/bbf_manouver.md)
- [`bff_eigenalaisy`](cartelle/bff_eigenalaisy.md)
- [`bff_longitudinal`](cartelle/bff_longitudinal.md)
- [`results`](cartelle/results.md)

Per aggiornare inventari e progressi dopo una modifica al progetto:

```bash
cd /home/nicomonzi/TESI
python3 report/aggiorna_report.py --note "descrizione sintetica del lavoro"
```

Lo script non modifica nulla fuori da `report`, non esegue simulazioni e non
cancella file. `.git` e `report` sono escluse dai conteggi per evitare che il
rapporto conti se stesso.
