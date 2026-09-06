# NASA MUTT — analisi modale, COUPMASS e convergenza MBDyn sotto 5 g

## Sintesi dei risultati

Il pacchetto contiene tre analisi collegate:

1. estrazione SOL 103 di 60 modi per il modal joint MBDyn;
2. sensibilità delle frequenze alla formulazione della massa tramite `COUPMASS`;
3. SOL 101 sotto gravità 5 g e sweep MBDyn da 1 a 54 modi elastici.

Risultati principali:

| Analisi | Risultato |
|---|---:|
| Modi Nastran estratti | 60 totali: 6 rigidi + 54 elastici |
| Intervallo elastico | 3.217–62.861 Hz |
| Minima distanza tra modi adiacenti | 0.0506 Hz |
| Coppie COUPMASS affidabili | 45 su 54 con MAC ≥ 0.90 |
| Scarto COUPMASS medio sulle coppie affidabili | 0.200% |
| Massimo scarto COUPMASS affidabile | 0.826% |
| Convergenza interna MBDyn, soglia 0.01 in | **30 modi elastici** |
| Errore MBDyn–Nastran minimo alle tip | 1.139 in |
| Errore MBDyn–Nastran con 54 modi | 1.170 in |
| Esito della soglia esterna a 1% | non raggiunta |

La conclusione centrale è che **30 modi elastici sono sufficienti per convergere
la soluzione MBDyn rispetto alla base completa a 54 modi**, ma il confronto con
Nastran non è ancora una validazione chiusa: l’errore esterno raggiunge un plateau
molto superiore alla soglia e quindi non dipende dal numero di modi.

## Struttura

```text
NASTRAN_SIMULATIONS/
├── 01_SOL103_60_MODES/
│   ├── MAIN/sol103_60_modes.bdf
│   ├── MAIN/mbdyn_modal_60.fem
│   ├── analyze_modal_spectrum.py
│   └── results/
├── 02_COUPMASS_STUDY/
│   ├── MAIN/sol103_coupmass_lumped.bdf
│   ├── MAIN/sol103_coupmass_coupled.bdf
│   ├── compare_coupmass.py
│   └── results/
└── 03_GRAVITY_5G/
    ├── nastran/MAIN/sol101_gravity_5g.bdf
    ├── mbdyn/mbdyn_modal.fem
    ├── mbdyn/gravity_5g_template.mbd
    ├── run_full_convergence.py
    ├── analyze_convergence.py
    └── plots/
```

Ogni caso Nastran possiede una copia autonoma della directory `BULK`.

## 1. Estrazione SOL 103 a 60 modi

Il deck `01_SOL103_60_MODES/MAIN/sol103_60_modes.bdf` usa:

```text
EIGRL ... 60
DISPLACEMENT(PLOT)=ALL
VECTOR(PLOT)=ALL
```

L’output completo su tutti gli 8527 nodi è necessario a `femgen`. Il file
generato e validato contiene:

- 8527 nodi FEM;
- 60 modi normali;
- matrice di rigidezza modale;
- diagonale della massa lumped.

I primi sei modi hanno frequenza quasi nulla e rappresentano i modi di corpo
rigido; il modal joint usa quindi i modi elastici Nastran 7–60.

![Spettro modale](01_SOL103_60_MODES/results/modal_spectrum_60.png)

La distanza minima di 0.0506 Hz mostra che lo spettro diventa molto denso nella
parte alta. In quella regione il solo numero del modo non è sufficiente per
seguirne l’identità: è necessario confrontare anche le forme modali.

## 2. Sensibilità a COUPMASS

Sono stati confrontati due SOL 103 identici salvo:

```text
PARAM,COUPMASS,-1   $ lumped mass
PARAM,COUPMASS,1    $ coupled/consistent mass
```

`COUPMASS` è un interruttore e non un coefficiente continuo. Un valore positivo
richiede la formulazione coupled, mentre il valore negativo predefinito usa la
formulazione lumped. Riferimento: [MSC Nastran documentation](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_SOL_400_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2022.1_SOL_400_Getting_Started_Guide.pdf).

Le forme sono abbinate sui 45 nodi del modal joint mediante MAC, imponendo anche
una finestra di frequenza del 15% per evitare associazioni tra bande remote.

- 45 coppie elastiche hanno `MAC ≥ 0.90`;
- su queste coppie lo scarto medio assoluto è 0.200%;
- lo scarto massimo affidabile è 0.826%, per il modo lumped 52 abbinato al
  coupled 50, circa 57 Hz;
- nove coppie sono a bassa confidenza, prevalentemente nella regione modale densa.

![Confronto COUPMASS](02_COUPMASS_STUDY/results/coupmass_comparison.png)

I primi modi elastici sono molto stabili. Ad esempio:

| Modo lumped | Frequenza [Hz] | Variazione coupled | MAC |
|---:|---:|---:|---:|
| 7 | 3.217 | 0.042% | 1.000 |
| 8 | 5.303 | 0.031% | 1.000 |
| 9 | 8.705 | 0.003% | 1.000 |
| 10 | 11.164 | 0.102% | 1.000 |
| 11 | 12.257 | 0.020% | 1.000 |

Lo script calcola anche copertura e partecipazione lungo lo span come indicatori
di modi globali. I candidati ad alta copertura includono i modi 8, 11, 20, 29,
35, 40–43, 48, 52 e 56. Questa classificazione è uno screening: i modi devono
essere confermati osservando le forme nell’OP2, soprattutto quando sono ravvicinati.

## 3. SOL 101 con gravità 5 g

Il modello statico usa:

- riferimento `990001`;
- RBE2 alla radice;
- SPC sulle componenti `123456` del riferimento;
- `GRAV = 1932 in/s²` nella direzione globale `-Z`;
- 45 nodi del modal joint stampati nel F06.

Gli spostamenti Nastran alle tip pre-winglet sono:

| Tip | Ux [in] | Uy [in] | Uz [in] |
|---|---:|---:|---:|
| Sinistra, 990020 | -0.3112 | -0.1136 | -13.8624 |
| Destra, 991020 | -0.3083 | +0.1126 | -13.7731 |

## 4. Sweep automatico MBDyn

L’intero workflow è contenuto in un solo comando:

```bash
python3 03_GRAVITY_5G/run_full_convergence.py --jobs 3
```

Lo script:

1. esegue `femgen` su OP2/MAT;
2. verifica 8527 nodi, 60 modi e i record di massa/rigidezza;
3. installa il FEM nel modal joint;
4. genera 15 casi con 1–54 modi elastici;
5. esegue MBDyn in parallelo con log separati;
6. legge il F06 SOL 101;
7. genera tabella, grafici e raccomandazione finale.

Lo sweep iniziale usa le basi:

```text
1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 30, 40, 50, 54
```

Per completare lo sweep esaustivo con ogni base intera da 1 a 54 si usa:

```bash
python3 03_GRAVITY_5G/run_remaining_54_bases.py
```

Lo script verifica i NetCDF, salta automaticamente i casi già conclusi e lancia
sequenzialmente soltanto quelli mancanti. Può essere interrotto e ripreso con lo
stesso comando. Per controllare lo stato senza generare o eseguire casi:

```bash
python3 03_GRAVITY_5G/run_remaining_54_bases.py --dry-run
```

Lo sweep esaustivo è stato completato e il CSV contiene ora tutte le basi intere
da 1 a 54. Al termine di eventuali nuovi run lo script ricrea automaticamente
CSV e grafico di convergenza.

La gravità è portata da 0 a 5 g in 2 s. È applicato il 20% di smorzamento modale
per estinguere rapidamente il transitorio senza modificare l’equilibrio statico.
La risposta è mediata sull’ultimo 20% della storia; il massimo residuo
picco-picco è circa 0.0044 in, inferiore alla soglia interna.

## 5. Convergenza interna MBDyn

Per isolare il solo errore di troncamento, ogni base è confrontata con quella più
ricca da 54 modi. La soglia interna scelta è 0.01 in sull’errore vettoriale massimo
alle due tip.

| Modi elastici | Errore vettoriale vs 54 modi [in] |
|---:|---:|
| 1 | 0.1701 |
| 10 | 0.0275 |
| 24 | 0.0124 |
| **25** | **0.0094** |
| 26 | 0.0103 |
| 27 | 0.0070 |
| 36 | 0.0102 |
| 37 | 0.0090 |
| 40 | 0.0050 |
| 50 | 0.0002 |
| 54 | riferimento |

![Convergenza interna](03_GRAVITY_5G/plots/gravity_5g_modal_truncation_error.png)

Versione vettoriale per LaTeX:
`03_GRAVITY_5G/plots/gravity_5g_modal_truncation_error.pdf`. Il grafico è
generato direttamente dal CSV mediante
`03_GRAVITY_5G/plot_modal_truncation_convergence.py`; il PNG è salvato a 450 dpi.

La prima base sotto 0.01 in è quella con **25 modi elastici**. La convergenza non
è strettamente monotona: 26 e 36 modi superano leggermente la soglia. La prima
sequenza di almeno tre basi sotto soglia inizia a 27 modi, mentre da **37 modi**
in poi tutti i casi più ricchi analizzati rimangono sotto 0.01 in.

## 6. Validazione esterna rispetto a Nastran

La soglia esterna è:

```text
max(0.01 in, 1% della massima deformazione tip Nastran) = 0.1387 in
```

Nessuna base la soddisfa. Con 54 modi:

| Tip | Errore vettoriale [in] | Errore relativo vettoriale | Errore Uz [in] | Errore relativo Uz |
|---|---:|---:|---:|---:|
| Sinistra | 1.1692 | 8.43% | 0.4365 | 3.15% |
| Destra | 1.1703 | 8.49% | 0.4466 | 3.24% |

![Errore rispetto a Nastran](03_GRAVITY_5G/plots/gravity_5g_nastran_validation_error.png)

### Errore assoluto lungo il semispan

Il confronto spanwise usa il semispan destro dal nodo di radice 990001 alla tip
991020 (`b = 0–168 in`), senza includere i nodi della winglet. I grafici mostrano
soltanto le basi con 1, 8, 25 e 54 modi elastici.

![Errore assoluto Ux](03_GRAVITY_5G/plots/gravity_5g_semispan_absolute_error_Ux.png)

![Errore assoluto Uy](03_GRAVITY_5G/plots/gravity_5g_semispan_absolute_error_Uy.png)

![Errore assoluto Uz](03_GRAVITY_5G/plots/gravity_5g_semispan_absolute_error_Uz.png)

Sono riportati anche gli errori assoluti delle tre componenti di rotazione,
confrontando `Phi` MBDyn con `Rx`, `Ry`, `Rz` del displacement vector Nastran;
tutti i valori angolari sono espressi in radianti.

![Errore assoluto Rx](03_GRAVITY_5G/plots/gravity_5g_semispan_absolute_error_Rx.png)

![Errore assoluto Ry](03_GRAVITY_5G/plots/gravity_5g_semispan_absolute_error_Ry.png)

![Errore assoluto Rz](03_GRAVITY_5G/plots/gravity_5g_semispan_absolute_error_Rz.png)

### Errori locali lungo il semispan

L'analisi corretta è in `03_GRAVITY_5G/local_spanwise_error`. Per evitare che
l'errore di una stazione contenga la deformazione accumulata nelle stazioni
precedenti, le traslazioni sono confrontate tramite la seconda derivata
spaziale `d²u/db²`, mentre le rotazioni tramite il gradiente locale `dtheta/db`.
Sono usate 18 stazioni interne per le curvature traslazionali e 19 segmenti per
i gradienti di rotazione, sempre senza winglet.

I grafici e il CSV si rigenerano con:

```bash
python3 03_GRAVITY_5G/local_spanwise_error/plot_local_spanwise_errors.py
```

### Indicatore integrale dell'errore spanwise

La cartella `03_GRAVITY_5G/spanwise_area_error` contiene l'indicatore
adimensionale ottenuto dall'area sotto la curva dell'errore composto. A ogni
nodo, l'errore angolare è convertito in spostamento equivalente tramite il
braccio locale `b`. Questo contributo è combinato con l'errore di traslazione;
l'area risultante è divisa per il quadrato della lunghezza del semispan ed
espressa come percentuale di `L = 168 in`. Il CSV riporta anche le aree separate
di Ux, Uy, Uz e delle rotazioni equivalenti Rx, Ry, Rz.

```bash
python3 03_GRAVITY_5G/spanwise_area_error/compute_spanwise_area_error.py
```

Il plateau dimostra che lo scarto non è dovuto alla quantità di modi. La causa
di coerenza più evidente è che il FEM modale SOL 103 è stato estratto usando il
file originale con **RBE3 alla radice**, mentre il SOL 101 di riferimento usa
**RBE2**. Inoltre la deformazione di circa 14 in è abbastanza grande da rendere
sensibile il confronto tra la cinematica MBDyn e una SOL 101 lineare.

Di conseguenza:

- **30 modi** sono sufficienti per la convergenza interna MBDyn;
- non è corretto dichiarare il modello validato a 1% contro questo SOL 101;
- per chiudere la validazione occorre estrarre un nuovo FEM modale con lo stesso
  RBE2 usato dal SOL 101 e ripetere lo sweep;
- se lo scarto permane, va confrontata la risposta MBDyn con una soluzione
  Nastran geometricamente non lineare oppure va ridotto il livello di g per una
  verifica strettamente lineare.

## Output numerici

- frequenze: `01_SOL103_60_MODES/results/modal_frequencies_60.csv`;
- matching COUPMASS: `02_COUPMASS_STUDY/results/coupmass_mode_matching.csv`;
- convergenza 5 g: `03_GRAVITY_5G/plots/gravity_5g_convergence.csv`;
- errori assoluti spanwise: `03_GRAVITY_5G/plots/gravity_5g_semispan_absolute_errors.csv`;
- log MBDyn: `03_GRAVITY_5G/mbdyn/results/*.run.log`.

Tutti i grafici usano la palette di tesi richiesta, linee continue e marker
soltanto per i punti discreti dello sweep.
