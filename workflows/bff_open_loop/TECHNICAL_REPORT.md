# Rapporto tecnico della ristrutturazione BFF

## Problema osservato

La run precedente a 65 m/s non era una prova BFF valida. Prima dello
spegnimento SAS la quota oscillava di molti metri; al rilascio si avevano
`Vz=-2.96 m/s` e `q=-15.27 deg/s`. Quota e Vz saturavano a +/-2 gradi,
pitch e q a +/-6 gradi. MBDyn divergeva a 19.86 s, prima del riaggancio.

La causa era un doppio percorso di controllo: quota e Vz modificavano sia il
riferimento di pitch sia direttamente WF1--WF3. Filtri, ritardi e saturazioni
facevano contrastare i due anelli. Inoltre lo smorzatore modale chiuso pari a
60 e un rate loop q troppo debole amplificavano il transitorio accoppiato.

## Modifiche al controllo

1. Quota e Vz agiscono solo in modo simmetrico su WF1--WF3.
2. Pitch e q sono separati e agiscono sui body flap.
3. Eliminato l'integrale del pitch durante l'ingresso al punto.
4. Ripristinati `Kp_pitch=-0.60*SAS_Q_SCALE` e
   `Kp_q=-0.75*SAS_Q_SCALE`.
5. Filtro Vz portato da 0.35 a 1 Hz per ridurre il ritardo di fase.
6. Aggiunto filtro q a 1 Hz per non inseguire il bending con i body flap.
7. Conservato il controllo roll/roll-rate differenziale.
8. Smorzamento modale chiuso schedulato: 16 fino a 65 m/s, 32 a 67.5 m/s,
   50 a 70 m/s. Il gate lo rende esattamente nullo durante SAS-off.

Una sonda chiusa a 65 m/s ha mostrato, nell'ultimo secondo prima di 8 s,
escursione quota 8.5 mm, `|Vz|max=0.0104 m/s` e
`|q|max=0.0366 deg/s`, contro la divergenza della configurazione precedente.

## Protocollo e sicurezza numerica

- Settling finale fino a 10.5 s.
- SAS-off limitato a 2.05 s.
- Hold simultaneo di tutte le dieci superfici.
- Solo il rap noto WF4 viene sommato al valore congelato.
- Riaggancio automatico e simulazione fino a 15.5 s.
- `run_case.py` controlla Vz, q, p e le deviazioni standard dell'ultimo
  secondo; un rilascio non stazionario arresta lo sweep.
- `analyse_open_loop.py` rifiuta i NetCDF interrotti prima del riaggancio.

Il tentativo a 70 m/s con smorzamento chiuso 16 e' stato respinto dal gate.
Anche i tentativi con 28 e 42 hanno mostrato assestamento insufficiente. Con
50 il rilascio finale e' `Vz=+0.00588 m/s`, `q=+0.00066 deg/s` e
`p=+0.00186 deg/s`.

## Aerodinamica e NASTRAN

- Tutti i 58 elementi usano Wagner/Theodorsen e collocazione a 3/4 di corda.
- La C81 efficace mantiene il drag NACA 0012 e corregge le pendenze globali:
  `CZ_alpha=0.108331/deg` e `Cm_alpha=-0.007550/deg`, contro NASTRAN
  `0.108327/deg` e `-0.007563/deg`.
- La verifica DUST FINE/VLM forniva `CZ_alpha=0.102012/deg` e
  `Cm_alpha=-0.009081/deg`; resta una verifica di tendenza, non la sorgente
  della frontiera flutter.
- La correzione ROM sul modo FEM 7 compensa i termini DLM/RFA 3-D mancanti.
  E' sempre presente in open loop e non e' una forza SAS.
- Il risultato SOL 145 raffinato attraversa zero a 65.7637 m/s TAS,
  equivalente a 60.5852 m/s EAS. La colonna VELOCITY dell'F06 e' KEAS e non
  deve essere confrontata direttamente con V_INF TAS.

## Identificazione e risultati

Il polo BFF e' selezionato dalla coerenza fra deformazione simmetrica delle tip
nel frame corpo, coordinata SWB1 e velocita' SWB1. q, alpha e pitch sono usati
per il candidato short-period. A 65 m/s il polo rigido coincide con quello
strutturale ed e' classificato come coupled BFF/short-period.

Risultati di controllo della ristrutturazione:

| TAS | f MBDyn | sigma MBDyn | f NASTRAN | sigma NASTRAN |
|---:|---:|---:|---:|---:|
| 57.5 | 2.3659 Hz | -3.1544 1/s | 2.3788 Hz | -3.1373 1/s |
| 65.0 | 1.9980 Hz | -0.4435 1/s | 2.0552 Hz | -0.2057 1/s |
| 70.0 | 1.6852 Hz | +1.2142 1/s | 2.0644 Hz | +0.9422 1/s |

Il segno cambia fra 65 e 70 m/s. Una interpolazione lineare preliminare da'
66.34 m/s TAS; lo sweep deve raffinare questo valore. La frequenza a 70 m/s e'
sottostimata: la ROM monomodale riproduce il cambio di stabilita', ma non e'
una sostituzione quantitativa completa del modello DLM/RFA multimodale.

## File operativi

- `run_case.py`: rendering, scheduling, run singola e gate sul trim.
- `analyse_open_loop.py`: audit, identificazione e grafici.
- `run_sweep.py`: griglia, tracking, raffinamento del crossing e riepilogo.
- `run_sweep.sh`: launcher WSL/Linux.
- `run_sweep_windows.bat`: launcher da Windows tramite WSL.
