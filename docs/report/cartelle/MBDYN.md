# Rapporto: `MBDYN`

## Ruolo e stato

È l'area storica di sviluppo MBDyn: 217 file e circa 784 MB. Contiene modelli,
esperimenti di controllo, output NetCDF e copie FEM. È **un archivio misto**:
alcune prove restano utili come provenienza, ma non è più il punto di ingresso
dei workflow correnti.

## Sottocartelle

- `MBDYN_MASS_CHECK`: verifiche di massa e FEM;
- `MBDYN_OK`: modelli strutturali/aerodinamici iniziali e script di animazione,
  frequenze e geometria;
- `MBDYN_SOL101`: confronto statico, numerosi NetCDF equivalenti e un grande
  testo risultati; è la porzione più pesante, circa 219 MB;
- `AUTOPILOT`, `PID`, `PID_test`, `PID_pitch`, `PID_5deg_35ms`: genealogia dei
  controllori. Gli ultimi due hanno runner e dashboard riproducibili per test
  pitch/roll a 35 m/s;
- `BBF`: modello BFF intermedio;
- `V-g_clamp`: sweep free-response e costruzione V-g/V-f mediante
  `analysiss.py`/`vgrun.py`.

## Codice

Gli script più vecchi sono procedurali e leggono direttamente MOV/NetCDF per
grafici e animazioni. I rami PID più recenti hanno entry point puliti,
docstring e dashboard. `V-g_clamp/analysiss.py` applica band-pass, stima
smorzamento/frequenza e interpola crossing. Le configurazioni MBDyn contengono
le prime versioni dei joint, aerobody, notch, PID e attuatori poi consolidati
nei workflow BFF.

## Duplicati e obsolescenza

La cartella ha 46 copie extra interne e un massimo teorico di 351 MB, quasi
tutto spiegato da FEM/NetCDF identici. I tre `xxx.nc`, `x.nc`, `symmetrical.nc`
e simili vanno confrontati col log/input prima dell'archiviazione.

Come linea di esecuzione, `AUTOPILOT`, `PID*`, `BBF` e `V-g_clamp` sono
**superseded** da `TEST`, `TRIM`, `bff_*` e infine `BFF_open_loop`. Conservare
README e casi minimi, migrare gli output unici in archivio e rimuovere cache;
non eliminare i FEM finché non esiste un registro canonico.

