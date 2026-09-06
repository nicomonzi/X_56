# Rapporto: `TEST`

## Ruolo e stato

È l'area di prototipazione che ha preceduto i workflow finali. Contiene 846 file
e circa 1,8 GB, la cartella di lavoro più grande escluso `.git`. È
**superseded come linea attiva**, ma contiene evidenze e algoritmi che spiegano
come si è arrivati ai modelli correnti.

## Sottocartelle

- `TRIM_PID_VELOCITY` (~1,7 GB): sweep 20–45 m/s con 26 NetCDF da circa 20 MB,
  file `.aer`, history CSV e analisi. `trim_velocity_sweep.py` prepara/esegue e
  `analyze_trim_sweep.py` costruisce riepiloghi e grafici;
- `BBF_AUTO_TRIM` (~27 MB): sweep 23–36 m/s, boundary/refinement, Hilbert,
  root-locus multimodale e dashboard. È il precursore dell'identificazione più
  rigorosa;
- `BBF_DUST_COUPLED` (~73 MB): adapter preCICE, runner coupling, sweep e
  analisi; contiene una build DUST locale da 33 MB e sei tentativi incompleti;
- `bbf` (~40 MB): singolo modello NASA-style con audit e post-processing.
- `FILTERS_AND_ACTUATORS.md`: riferimento tecnico su notch, low-pass, Tustin e
  taratura dei guadagni.

## Valore tecnico

Questa cartella documenta l'evoluzione da PID e trim automatico verso
identificazione BFF, gestione delle saturazioni e coupling DUST. Gli script
contengono molte funzioni ancora valide, ma versioni più recenti hanno corretto
la separazione dei loop, lo stato al rilascio, la trasparenza del SAS e la
robustezza dell'identificazione.

## Obsolescenza e pulizia

`build_dust`, `__pycache__` e `incomplete_attempt_*` sono rigenerabili. Gli
output dello sweep trim sono archiviabili dopo aver conservato summary, input,
versione solver e almeno i casi rappresentativi. Non cancellare in blocco:
alcuni NetCDF possono essere l'unica prova delle campagne storiche. Spostare
l'intera area sotto `archive/prototypes_test` con un indice che mappi ogni ramo
al successore: `TRIM`, `BFF_open_loop` o `BFF_DUST_55`.

