# X-56A BFF — bounded full-surface-hold SAS-off test

Questa cartella esegue una prova locale simile al protocollo NASA usato per
avvicinarsi al body-freedom flutter: ingresso stabilizzato nel punto di prova,
breve rap simmetrico, tutte le superfici congelate, osservazione libera e
riattivazione automatica del SAS.

## Architettura finale

- Aerodinamica sezionale: Wagner/Theodorsen su tutti i 58 elementi, polar
  `INCLUDE/x56_effective.c81` e collocazione a 3/4 di corda.
- La polar statica C81 riproduce le pendenze globali NASTRAN di portanza e
  momento. Il drag NACA 0012 viene conservato.
- Una correzione ROM sul solo modo FEM 7 rappresenta i termini DLM/RFA 3-D
  mancanti nel modello sezionale. Questa forza aerodinamica rimane attiva
  durante SAS-off e non e' un guadagno di controllo.
- Il controllo di traiettoria usa quota e velocita' verticale soltanto sulle
  superfici alari simmetriche WF1--WF3.
- Pitch e pitch-rate agiscono separatamente sui body flap. `q` e' filtrato a
  1 Hz per non inseguire direttamente il bending.
- Roll angle e roll-rate sono controllati differenzialmente su WF1--WF2.
- Sopra 65 m/s uno smorzatore modale schedulato rappresenta il sistema attivo
  di flutter suppression necessario per entrare e recuperare il punto. Il suo
  comando e' esattamente zero durante la finestra open-loop.

Quota/Vz non modificano piu' anche il riferimento di pitch: quel doppio
percorso generava il ciclo limite rigido osservato nella vecchia run a 65 m/s.

## Sequenza

1. SAS chiuso e assestamento fino a `10.5 s`;
2. verifica automatica dello stato al rilascio (`Vz`, `q`, `p` e dispersioni);
3. hold di BFL/BFR e WF1--WF4 al valore dell'ultimo campione chiuso;
4. piccolo rap simmetrico WF4 centrato sulla frequenza NASTRAN;
5. identificazione tip/SWB1/SWB1dot con tutti i feedback nulli;
6. riattivazione SAS dopo `2.05 s` e recupero fino a `15.5 s`.

Se il trim non soddisfa le soglie, `run_case.py` termina con errore e lo sweep
si arresta: un caso rilasciato con deriva non puo' diventare un risultato BFF.
Anche `analyse_open_loop.py` rifiuta esplicitamente un NetCDF incompleto.

## Identificazione

Il BFF viene selezionato come coppia complessa persistente e coerente fra:

- deformazione verticale simmetrica delle tip nel frame del velivolo;
- coordinata modale SWB1 (modo FEM 7);
- velocita' modale SWB1.

Traslazione e rotazione rigide vengono eliminate dalla deformazione delle tip.
`q`, angolo d'attacco effettivo e pitch identificano separatamente lo
short-period e indicano se questo partecipa allo stesso polo del bending.

## Verifiche eseguite

I punti 57.5 e 65 m/s sono stati verificati con rilascio a 10.0 s durante il
tuning; il punto 70 m/s usa il valore finale di 10.5 s. Lo sweep rigenera tutti
i punti con la configurazione finale e sostituisce questi numeri preliminari.

| TAS [m/s] | Stato rilascio | BFF MBDyn | SOL 145 punto 7 | Esito |
|---:|:---|:---|:---|:---|
| 57.5 | `Vz=+0.0161 m/s`, `q=-0.0150 deg/s` | `f=2.3659 Hz`, `sigma=-3.1544 1/s` | `2.3788 Hz`, `-3.1373 1/s` | stabile, ottimo accordo |
| 65.0 | `Vz=-0.00536 m/s`, `q=-0.0119 deg/s` | `f=1.9980 Hz`, `sigma=-0.4435 1/s` | `2.0552 Hz`, `-0.2057 1/s` | stabile vicino al confine |
| 70.0 | `Vz=+0.00588 m/s`, `q=+0.00066 deg/s` | `f=1.6852 Hz`, `sigma=+1.2142 1/s` | `2.0644 Hz`, `+0.9422 1/s` | instabile; frequenza ROM sottostimata |

Il cambio di segno di `sigma` identifica un confine fra 65 e 70 m/s, coerente
con SOL 145 (`65.764 m/s TAS`). L'interpolazione lineare dei soli punti 65 e
70 fornisce un valore preliminare MBDyn di circa `66.34 m/s TAS`; e' lo sweep
che deve calcolare il valore definitivo. Il punto alto dimostra il limite della
correzione monomodale: riproduce il passaggio stabile/instabile ma non ancora
la frequenza DLM quantitativa a 70 m/s. Lo sweep raffina il crossing MBDyn;
non deve essere presentato come una nuova soluzione DLM indipendente.

## Sweep 50--70 m/s

Da WSL/Linux:

```bash
cd /home/nicomonzi/X_56/workflows/bff_open_loop
./run_sweep.sh
```

Oppure da Windows avviare `run_sweep_windows.bat`. Il launcher esegue:

```bash
python3 run_sweep.py --start 50 --stop 70 --step 2.5 \
  --refine-tolerance 0.25 \
  --output /mnt/c/Users/Utente/Desktop/BFF_open_loop \
  --clean --overwrite
```

`--clean` elimina soltanto file prodotti da precedenti sweep, non altri file
dell'utente. I risultati vengono salvati in
`C:\Users\Utente\Desktop\BFF_open_loop`. Lo sweep non viene avviato
automaticamente.

## Riferimenti NASTRAN

Il parser usa prima il risultato raffinato:
`/home/nicomonzi/ZENO/X56_NASTRAN_BFF_REFINEMENT/FLUTTER_TEST/x56_bff_refined.f06`.
La colonna `VELOCITY` dell'F06 e' KEAS; viene convertita in TAS usando `VREF`.
Il crossing del punto 7 e' `65.7637 m/s TAS`, pari a `60.5852 m/s EAS`.
