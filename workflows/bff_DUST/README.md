# BFF open-loop con DUST: singola velocita'

Workflow separato da `bff_open_loop` e `coupled_dust`: non modifica i casi
originali e non esegue sweep. Gli ingressi vengono generati in una nuova
sottocartella `runs/` per ogni comando; una cartella esistente non viene sovrascritta.

## Configurazione

- FEM60: 25 modi elastici, indici 7--31, 3.217134--36.970058 Hz. Modi rigidi FEM
  esclusi; moto rigido gestito dal nodo modale MBDyn. Smorzamento strutturale 1%.
- Mesh FINE selezionata: `nelem_chord=30`, span `[2,3,4,3,9,40]`,
  7320 pannelli, 7503 vertici, dieci superfici mobili.
- `dt=0.002 s`: 13.525 campioni/periodo del modo piu' alto; questa scelta non
  sostituisce una verifica di convergenza temporale del flutter accoppiato.
- IPS coerente: lunghezze in, forze lbf, momenti lbf in, tempo s,
  densita' 9.7284e-8 lbf s^2/in^4. Velocita' nel JSON in m/s, convertita in in/s.
- Unico punto iniziale: **66.34 m/s TAS** = 2611.811024 in/s. E' il crossing
  MBDyn *preliminare* documentato nel README di `bff_open_loop`, non una velocita'
  di flutter DUST gia' verificata. Non e' disponibile qui il crossing definitivo
  dello sweep. Il parametro e la provenienza sono espliciti in `case.json`.

## Controlli e sequenza

### Scia limitata (aggiornamento successivo al primo smoke)

`n_wake_panels=2`: due file di pannelli di scia prima della conversione in
vortoni. I vortoni vengono eliminati oltre `x=581 in` mediante
`particles_box_max`: 15 corde da 24 in dopo `x=221 in`, arrotondamento del
bordo d'uscita piu' a valle della mesh (220.725 in). La corda di riferimento
24 in e' quella di `bff_open_loop/dust_static_validation.py`, non la corda
del centrocorpo. Sono 360 in = 9.144 m di scia dal riferimento scelto.
Il piano e' fisso globale: non segue il velivolo e non e' un limite di eta'
dei vortoni. Restano i precedenti limiti laterali/verticali +/-1500 in e
il limite a monte x=-1000 in. La capacita' e' 50000 vortoni, indipendente
dalla durata; questa capacita' non sostituisce la cancellazione spaziale.

Il tempo misurato di 365.119 s e il primo smoke si riferiscono alla vecchia
scia di 40 file, senza questo taglio a 15 corde. Non sono un benchmark della
nuova configurazione. Il taglio puo' modificare i carichi e il flutter:
la sensibilita' alla lunghezza della scia resta da verificare.

### Controllori

Sono rimossi **soltanto YAW_PID e R_PID**, e il vincolo del nodo base blocca
la terza componente di rotazione relativa al riferimento del giunto.
Restano ALT, VZ, PITCH, Q, ROLL, P e VY e lo smorzatore modale del caso originale.
VY conserva il percorso originale: il limite direzionale e' zero, percio'
questo PID resta presente ma non genera una deflessione direzionale.

Filtri ALT/VZ/Q/VY e attuatori vengono ridiscretizzati a ogni generazione usando
`dt_s`, i cutoff e `actuator_tau_s` nel JSON. Nessun coefficiente a dt=0.01 s
viene riutilizzato. I guadagni continui PID non sono coefficienti di filtro:
restano **valori iniziali da tarare con DUST**, con moltiplicatori indipendenti
`controller_gain_scale`; lo smorzatore modale e' configurabile separatamente.
Non viene dichiarata una taratura DUST a partire da un test di pochi passi.

Il caso completo dura 15.5 s:

1. Assestamento con controllori fino a 10.5 s; controllo dello stato prima del rilascio.
2. Hold delle dieci superfici e smorzatore modale nullo.
3. Rap simmetrico WF4 da 10.55 s, ampiezza 0.2 gradi, durata 0.742/2.057 s.
4. Osservazione libera dopo la coda del filtro (almeno 8 costanti di tempo).
5. Recupero con riattivazione SAS a 12.55 s.

Non e' una ricerca di trim e non include una manovra di pull-up. L'assestamento
serve a rendere interpretabile il rilascio open-loop. Se Vz, q, p o le loro
dispersioni superano le soglie, il launcher interrompe il caso prima del rilascio
e salva `release_check.json`: occorre tarare i controllori non-yaw.
Le superfici vengono congelate durante la finestra aperta, anche se gli stati
interni dei PID continuano a evolvere, come nel workflow originale.

Sono eliminati aerodinamica sezionale MBDyn, polar C81 e correzione ROM DLM.
Resta la gravita'; gli unici carichi aerodinamici arrivano da DUST via preCICE.
Lo smorzatore modale rimasto e' una forza di controllo, non aerodinamica.

## Accoppiamento e riproducibilita'

preCICE serial-implicit scambia posizione, rotazione, velocita', velocita'
angolare, forze e momenti su 59 nodi. I nodi ausiliari di cerniere adiacenti
che prima coincidevano sono arretrati di 0.01 in lungo la rispettiva cerniera,
sia in MBDyn sia in DUST: elimina i pareggi nearest-neighbor senza cambiare
i vertici della superficie aerodinamica. I primi 39 nodi rimangono invariati.

Ogni run salva ingressi, adapter, configurazione effettiva, hash del FEM e
audit dei modi/filtri. FEM e profili condivisi restano riferimenti assoluti
al repository: non spostare il caso senza rigenerarlo.

Il test tecnico verifica lettura del FEM, XML, geometria, iterazioni implicite,
valori finiti, vincolo yaw e hold. Non certifica il flutter, la convergenza
temporale o l'accuratezza fisica. L'adapter e il binario DUST sono quelli del
workflow `coupled_dust`; in quell'accoppiamento e' stato osservato uno scarto
di un dt tra etichette temporali delle forze H5 e CSV. Prima di usare i poli
per un confronto quantitativo di flutter occorre verificare anche la fase
dei carichi nell'ordine temporale dello scambio; la sola convergenza implicita
non dimostra questa proprieta'.

## Comandi

Da `/home/monzani/X_56/workflows/bff_DUST`:

```bash
python3 -m unittest discover -s tests -v
python3 run_case.py --prepare
python3 run_case.py --check
python3 run_case.py --smoke --threads 12
python3 run_case.py --run --threads 12
python3 run_case.py --analyse runs/NOME_RUN
```

`--smoke` comprime tutte le fasi in 0.04 s con la stessa mesh e gli stessi
25 modi: non verifica assestamento o flutter. `--run` avvia il caso lungo:
il costo e' elevato e aumenta con la scia. Non viene avviato automaticamente.
`--config FILE` permette una taratura esplicita mantenendo il caso a singola
velocita'; `--output CARTELLA_NUOVA` sceglie la destinazione.

Output: `analysis.json`, `response.csv`, `coupled_response.csv`, NetCDF MBDyn,
HDF5 DUST e serie `paraview/x56.pvd`, `x56_wpan.pvd`, `x56_wpart.pvd`.
La GUI ParaView richiede una sessione grafica disponibile.

I binari predefiniti sono MBDyn `/usr/local/mbdyn/bin/mbdyn` e DUST
`/home/monzani/dust-patched/build-user/bin/`. Per cambiarli usare le variabili
`MBDYN_BIN`, `DUST_BIN`, `DUST_PRE_BIN`, `DUST_POST_BIN`, `MBDYN_PYTHON_PATH`
oppure un file locale `machine.env` con righe `NOME=/percorso`.
Il launcher termina solo i propri processi in caso di errore/interruzione
e interrompe uno stallo senza avanzamenti dopo 900 s (configurabile).
