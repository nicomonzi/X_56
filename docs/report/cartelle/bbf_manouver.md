# Rapporto: `bbf_manouver`

## Ruolo e stato

È il primo workflow completo per BFF durante manovre level, pull-up e roll.
Contiene 44 file e circa 245 MB; 208 MB sono output NetCDF/log. È un
**precursore superseded** da `BFF_maneuver_envelope` e
`MANOUVER_STIFNESS`, ma quattro configurazioni a 60,8421 m/s risultano
documentate come validate.

## Codice e modello

- `run_manouver.py`: seleziona modalità, burst e velocità e genera un input
  descrittivo senza sovrascrivere il base;
- `run_sweep.py`: varia soltanto V_INF e conserva output per stem;
- `analyse_bff.py`: Matrix Pencil, ricostruzione dei comandi, residui
  aerodinamici e grafici; `analyse_bff copy.py` è una copia quasi/totalmente
  ridondante condivisa con altri rami;
- `INCLUDE/`: FEM, nodi, joint, aerobody e polar;
- `output/`: cinque NetCDF fra 25,6 e 47,7 MB più log e input renderizzati.

Il modello acquisisce il trim per 12 s, rilascia il volo, entra in manovra fra
14 e 22 s e osserva il recupero fino a 28 s. X/Y restano vincolati; Z, roll,
pitch e yaw sono liberi. Il SAS rimane sempre attivo, perciò il caso non misura
un polo open-loop come i workflow successivi.

## Obsolescenza e azioni

Archiviare come baseline storica, mantenendo i quattro casi validati e i loro
input. Deduplicare lo script `copy.py` dopo confronto col file corrente.
Migrare le conclusioni ancora citate in un indice di provenienza, poi evitare
nuove esecuzioni da questa cartella.

