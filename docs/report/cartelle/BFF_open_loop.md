# Rapporto: `BFF_open_loop`

## Ruolo e stato

È la baseline MBDyn steady più matura per l'identificazione BFF con rilascio
controllato del SAS. Contiene 35 file, circa 56 MB, tutti tracciati. È
**canonica** e viene importata dai workflow di manovra; non è obsoleta.

## Modello

Il modello usa 58 elementi aerodinamici sezionali Wagner/Theodorsen, una C81
effettiva calibrata sulle pendenze globali Nastran e una correzione ROM DLM sul
modo FEM 7. I modi rigidi 1–6 del FEM non sono usati nel modal joint; i gradi di
corpo sono rappresentati dal floating frame. Il controllo separa quota/Vz,
pitch/q, roll/p e smorzamento modale schedulato. Durante la finestra open-loop i
feedback e il damper sono realmente nulli e tutte le dieci superfici sono
congelate, salvo il RAP noto.

## Codice

- `run_case.py`: render di un caso, scheduling DLM/SAS, validazione dello stato
  al rilascio ed esecuzione MBDyn;
- `run_sweep.py` e `run_sweep.sh`: griglia 50–70 m/s e raffinamento del crossing;
- `analyse_open_loop.py`: oltre mille righe di controlli, ricostruzione comandi,
  separazione moto rigido/flessibile, tabelle e grafici;
- `modal_identification.py`: Matrix Pencil multimodale, clustering e tracking
  dei candidati BFF/short-period;
- `run_diagnostics.py`: otto casi diagnostici e continuità del ramo;
- `nastran_flutter_reference.py`: parser SOL 145, conversione KEAS→TAS e
  interpolazione del crossing;
- `build_x56_c81.py`, `aero_static_validation.py` e
  `dust_static_validation.py`: calibrazione/validazioni aerodinamiche.

Gli input MBDyn sono divisi fra file principali e `INCLUDE/` per nodi, joint,
forze, aerobody, FEM e polar. `TECHNICAL_REPORT.md` registra il difetto del
doppio percorso quota–pitch e la sua correzione; `TUNING_REPORT.md` documenta
la taratura.

## Risultati

I punti documentati sono: 57,5 m/s stabile (`f=2,366 Hz`, `sigma=-3,154 1/s`),
65 m/s stabile vicino al confine (`f=1,998 Hz`, `sigma=-0,444 1/s`) e 70 m/s
instabile (`f=1,685 Hz`, `sigma=+1,214 1/s`). L'interpolazione preliminare 65–70
fornisce circa 66,34 m/s, rispetto a 65,764 m/s TAS SOL 145.

La frequenza a 70 m/s è sottostimata; il modello riproduce il cambio di
stabilità ma non costituisce una nuova DLM indipendente. Lo sweep definitivo è
manuale e scrive fuori dal repository.

## Azioni

Mantenere questa cartella come workflow di riferimento. Centralizzare il FEM da
57,4 MB senza rompere gli include, salvare manifest dei risultati esterni e
aggiungere test sintetici unitari per `modal_identification.py`. Le copie C81 e
FEM non vanno eliminate prima della migrazione dei path.

