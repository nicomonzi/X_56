# Audit dei percorsi di comando

## Modello precedente `bff_longitudinal`

Ogni percorso capace di modificare una superficie e' stato tracciato:

| Percorso | Sensore / origine | Superficie nel modello precedente | Stato durante la vecchia identificazione |
|---|---|---|---|
| trim | fit dipendente da densita' e `V_INF` | tutte; bias distinto sui body flap | attivo e costante |
| quota | `Z_CG` -> `ALT_PID` | riferimento del `PITCH_PID`; ramo lift predisposto ma moltiplicato per zero | attivo |
| velocita' verticale | `XP[3]` -> `VZ_PID` | riferimento del `PITCH_PID`; ramo lift predisposto ma moltiplicato per zero | attivo |
| pitch | `E[2]` -> `PITCH_PID` | BFL/BFR simmetrici | attivo, attraverso quattro band-stop |
| pitch rate | `Omega[2]` -> `Q_PID` | BFL/BFR simmetrici | attivo, notch sensore + quattro band-stop |
| roll | `E[1]` -> `ROLL_PID` | WF1/WF2 differenziali | attivo |
| roll rate | `Omega[1]` -> `P_PID` | WF1/WF2 differenziali | attivo con notch |
| yaw | `E[3]` -> `YAW_PID` | BFL/BFR differenziali | attivo |
| yaw rate | `Omega[3]` -> `R_PID` | BFL/BFR differenziali | attivo con notch |
| velocita' laterale | `XP[2]` -> `VY_PID` | BFL/BFR differenziali | attivo |
| safety | deadband su pitch e q | BFL/BFR simmetrici | armato quattro secondi dopo il burst |
| burst | sinusoide Hann a frequenza fissata | WF4L/WF4R simmetrici | attivo durante l'eccitazione |
| attuatori | Tustin primo ordine, saturazione correzione | longitudinal, lateral, directional, lift | sempre attivi |
| sample-and-hold | assente | nessuna | assente |

Quindi nella vecchia finestra non esisteva true open loop: quota, Vz, pitch, q,
roll, p, yaw, r e Vy potevano continuare a cambiare le superfici; il safety
poteva intervenire piu' tardi. I notch riducevano ma non eliminavano il
feedback nella banda BFF.

## Caso attuale `BFF_open_loop`

- Il trim resta un bias fisico separato.
- Durante settling e recupero quota/Vz agiscono sia sul riferimento pitch sia, con
  banda bassa, simmetricamente su WF1--WF3. Pitch/q agiscono sulle body flap.
  Roll-angle e roll-rate agiscono differenzialmente su WF1/WF2; il
  direzionale resta disconnesso.
- Non esistono notch, safety controller, manovra, thrust o airspeed hold.
- A `SAS_OFF_START=7.00 s` i sample-and-hold congelano senza salto tutte le
  dieci superfici: BFL/BFR e WF1--WF4. Nessun feedback raggiunge il plant.
  WF4 applica il solo input noto, un doublet simmetrico da 0.20 deg attraverso
  l'attuatore Tustin (`tau=0.01 s`, `dt=0.01 s`), sommato al valore congelato.
- Da `IDENTIFICATION_START=7.25 s` a `IDENTIFICATION_END=8.00 s` il rap e'
  terminato e il canale BFF resta aperto. L'audit richiede costanti tutte le
  superfici, controlla la continuita' al rilascio e impone nulli tutti i
  feedback applicati e il loro leakage alla frequenza candidata.
- A `SAS_ON_START=8.05 s` i sample-and-hold tornano a seguire i comandi live e
  il SAS recupera automaticamente quota e assetto. I PID grezzi restano sempre
  registrati per audit.
- L'audit open-loop usa i comandi effettivi ricostruiti: tutte le superfici
  devono restare entro `0.001 deg`. La rotazione relativa dei nodi e' conservata
  come diagnostica separata, perche' contiene anche la deformazione elastica e
  non deve essere scambiata per feedback residuo.
- L'identificazione non usa `q`: la deformazione simmetrica delle tip viene
  calcolata nel frame del floating body, rimuovendo heave e rotazione rigida,
  e deve condividere il polo con SWB1 e SWB1dot.
- Tutti i 58 elementi usano `x56_effective.c81`; le pendenze globali sono
  verificate contro NASTRAN dopo il trasferimento del momento all'effettivo CG.

## Base modale

I modi rigidi FEM 1--6 sono esclusi perche' rappresentati dal floating frame.
Con `dt=0.01 s` restano attivi i modi elastici FEM 7--12: 3.2171, 5.3027,
8.7051, 11.1640, 12.2571 e 12.7589 Hz. Il postprocessore legge direttamente
coordinate e forme modali dal FEM, confronta le componenti verticali delle due
ali e identifica SWB1 come il modo simmetrico a frequenza piu' bassa; nel FEM
attuale risulta FEM 7, con residuale di simmetria 0.00369.
