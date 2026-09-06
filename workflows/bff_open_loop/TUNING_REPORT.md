# Final 62.5 m/s bounded SAS-off report

## Configurazione finale

- SAS aperto da 7.00 a 8.05 s con sample-and-hold comune a BFL/BFR e WF1--WF4.
- I controllori quota/Vz, pitch/q e roll/p continuano a evolvere internamente,
  ma nessuna loro uscita raggiunge le superfici durante SAS-off.
- Rap simmetrico WF4: 7.05--7.20 s, ampiezza 0.20 deg.
- Identificazione BFF: 7.25--8.00 s, usando tip simmetrica nel body frame,
  SWB1 e SWB1dot. Identificazione short-period separata con `q`, angolo
  d'attacco effettivo e pitch; questi canali non selezionano il ramo BFF.
- C81 effettivo X-56 calibrato sulle pendenze full-aircraft NASTRAN al CG.
- Riferimento di pitch aeroelastico schedulato sulla velocita' tra 50 e 70 m/s;
  superfici e coordinate modali restano sul feed-forward fisico originale.
- Roll PID: Kp=-0.80, Ki=-0.080; rate PID Kp=-0.40.
- Lift hold lento e filtrato: `LIFT_HOLD_GAIN=-0.58`.

## Risultato 62.5 m/s

- `OPEN_LOOP_VALID=True`, `ID_VALID=True`, trim stazionario.
- Rilascio: theta=0.87527 deg, q=-0.01042 deg/s, Vz=-0.000393 m/s,
  roll=-0.01690 deg, p=0.00107 deg/s.
- Quota prima del rap: -0.016 mm; alla riaccensione: -7.48 mm.
- Perdita massima durante SAS-off: 7.35 mm; massima complessiva: 10.91 mm.
- Roll massimo: 0.01690 deg; p massimo: 0.01279 deg/s.
- Stato finale: errore quota=-5.77 mm, Vz=-0.000073 m/s,
  roll=-0.000079 deg, p=0.000247 deg/s.
- BFL/BFR e WF1--WF3 sono costanti a precisione numerica. L'escursione WF4
  residua dopo il rap e' 0.000549 deg, sotto la tolleranza 0.001 deg.
- Tutti i feedback applicati e il leakage in banda BFF sono esattamente nulli.
- Polo BFF locale: f=2.68720 Hz, sigma=-3.54841 1/s, damping ratio=0.20567;
  correlazione tip--SWB1=-0.999982.
- Candidato q/alpha/pitch: circa 8.63 Hz, sovrapposto al modo FEM simmetrico 9
  secco (8.70512 Hz): non viene dichiarato short-period isolato.
- L'ampiezza armonica residua di p a circa 5.5 Hz e' 0.000205 deg/s, contro
  0.03478 deg/s con il precedente rate gain: riduzione superiore al 99%.

## Aerodinamica

Le pendenze rigide MBDyn/C81 (`CZ_alpha=0.108331/deg`,
`Cm_alpha=-0.007550/deg`) coincidono con NASTRAN trasferito al CG
(`0.108327/deg`, `-0.007563/deg`). La prova DUST FINE/VLM con 2142 pannelli e
0.4 s di sviluppo scia restituisce `0.102012/deg` e `-0.009081/deg`.
Il risultato DUST e' compatibile come controllo indipendente, ma non viene
usato per overfitting perche' mancano le winglet e l'audit della mesh riporta
14.59% di variazione MEDIUM--FINE sul momento.

## Stato degli output

Le run di tuning e validazione sono state scritte fuori dalla cartella del
modello. Lo sweep 50--70 m/s e' predisposto in `run_sweep.py`, salva per default
in `C:\Users\Utente\Desktop\BFF_open_loop` e non viene avviato automaticamente.
