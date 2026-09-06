# MANOUVER_STIFNESS

## Rettifica della revisione del 5 settembre 2026

Il precedente gate `physical_prestress_resolved_robust_to_dlm_distribution_build_rom`
non è validato: entrambe le SOL 103 a 1.6 g contengono una radice negativa
esclusa dal vecchio postprocess con la soglia di 0.1 Hz. Inoltre il recupero
del precarico usa coefficienti DLM diversi da quelli dei deck MBDyn eseguiti.
I risultati di integrazione restano disponibili, ma non dimostrano da soli
stabilità del precarico né robustezza alla distribuzione DLM.

La revisione completa di fisica, implementazione e campagna minima proposta
è in [REVIEW_MODAL_NASTRAN_DLM.md](REVIEW_MODAL_NASTRAN_DLM.md).
L'audit separato è riproducibile con `python3 review_existing_results.py`;
i precedenti JSON/CSV e i deck sono conservati, non corretti silenziosamente.

## Esito e scopo

Questa cartella contiene input e strumenti pronti per verificare, con un unico
metodo time-domain, se la richiamata modifica l'onset del body-freedom flutter
(BFF) e se il risultato è abbastanza sensibile alla rigidezza strutturale da
richiedere un vero modello prestressed.

La conclusione dell'audit è precisa: **il modal joint MBDyn attuale non include
lo stress-stiffening**. Include la massa distribuita e gli accoppiamenti
inerziali del corpo flessibile, ma usa una matrice elastica modale costante. Il
supporto MBDyn per la rigidezza geometrica esiste, però richiede le matrici di
carico unitario nel `RECORD GROUP 19` del file FEM e, nella versione installata,
il ramo con automatic differentiation. Il FEM X-56 contiene soltanto i gruppi
1–11, il modello non abilita l'automatic differentiation e la `femgen` locale
non genera il gruppo 19.

Per questo non ho costruito un falso prestress da 5 g. Ho preparato:

1. una campagna principale a rigidezza lineare costante;
2. una convergenza temporale coerente anche per filtri e attuatori;
3. uno screening parametrico della rigidezza del modo BFF, usato come gate
   quantitativo prima di spendere lavoro in un ROM stress-stiffened.

Le **28 traiettorie MBDyn ottimizzate** della campagna sono state eseguite e
analizzate. Sono state poi eseguite due sole shadow aggiuntive per recuperare i
carichi fisici a 1 g e 1.6 g; non è stato ripetuto lo sweep. Input, risultati,
carichi e deck Nastran sono salvati in
`C:\Users\Utente\Desktop\BFF_PULLUP_V2`. Nastran non è installato in questo
ambiente: le quattro SOL 103 frozen-time e la baseline supportata sono state
successivamente eseguite dall'utente su Zeno, con i limiti descritti nella rettifica.
La matrice estesa precedente è conservata in `runs_superseded_full_matrix` e
non va eseguita.

## Audit fisico del modal joint

Per un nodo FEM `P`, la cinematica modale floating-frame è, in forma
compatta,

\[
x_P=x_0+R_0\left(f_P+\Phi_Pq\right).
\]

Il corpo può quindi avere grandi traslazioni e rotazioni rigide, mentre la
deformazione elastica rappresentata da `q` resta piccola. Nel modello base
l'energia elastica linearizzata è

\[
\delta W\simeq \delta q^T K_{qq}q,
\]

con `K_qq` costante. La presenza del `RECORD GROUP 11` permette a MBDyn di
ricostruire proprietà di massa e invarianti: questo rende non banale la
dinamica inerziale del corpo flessibile, ma **non** trasforma `K_qq` in una
rigidezza dipendente dal carico.

Con il `RECORD GROUP 19`, MBDyn può invece usare

\[
\delta W\simeq\delta q^T\left(K_{qq}+K_{qq,geo}\right)q,
\]

dove `K_qq,geo` è una combinazione delle accelerazioni rigide, della
gravità, dei termini angolari e delle reazioni forza/momento ai nodi di
interfaccia. Le matrici richieste non sono le normali matrici `M_qq` e
`K_qq` esportate da `femgen`: sono matrici geometriche ottenute da campi di
sforzo per carichi unitari.

L'audit riproducibile è:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 audit_modal_joint.py
```

Il risultato completo è in `audit/modal_joint_audit.json`. Dati principali:

- 8527 nodi FEM e 60 modi mass-normalizzati;
- modi elastici selezionati: FEM 7–12;
- frequenze: 3.2171, 5.3027, 8.7051, 11.1640, 12.2571 e 12.7589 Hz;
- massa modale unitaria e matrici modali diagonali;
- gruppi presenti: 1–11;
- gruppo 19: assente;
- automatic differentiation nel modello corrente: assente;
- stato risultante: `linear_reduced_stiffness_only`.

I riferimenti locali verificati sono:

- `/home/nicomonzi/src/mbdyn/manual/tecman/tecman-modal.tex`, sezione
  `Quasi static corrections for stress stiffening`;
- `/home/nicomonzi/src/mbdyn/manual/input/modal-fem-format.tex`, definizione e
  tag del `RECORD GROUP 19`;
- `/home/nicomonzi/src/mbdyn/mbdyn/struct/modalad.cc`, applicazione dinamica
  delle matrici geometriche;
- `/home/nicomonzi/src/mbdyn/mbdyn/struct/modal.cc`, lettura del gruppo 19 e
  selezione del ramo `ModalAd`;
- `/home/nicomonzi/src/mbdyn/utils/femgen.f90`, che nella build locale scrive
  soltanto i gruppi fino all'11.

## Modello aeroelastico conservato

La campagna usa il modello validato in `BFF_open_loop` e il generatore di
manovra in `BFF_maneuver_envelope`. Restano invariati:

- i sei modi elastici FEM 7–12 e lo smorzamento modale strutturale 1%;
- la correzione calibrata DLM applicata al modo 7, fissata alla TAS nominale
  durante il rendering, non aggiornata con la TAS istantanea;
- il SAS e il damper modale durante ingresso e recupero;
- lo spegnimento del SAS per 2.05 s;
- il blocco di tutte le superfici al valore presente al rilascio;
- il rap simmetrico WF4 di 0.20° solo nella traiettoria excited;
- gravità reale pari a 1 g.

Il termine DLM non è un comando di controllo: rimane attivo durante SAS-off.
Il damper modale SAS è invece nullo nella finestra aperta.

Ogni condizione contiene due run con la stessa traiettoria:

- `shadow`: nessun rap BFF;
- `excited`: rap BFF noto.

Il segnale identificato è

\[
\Delta q_7=q_{7,excited}-q_{7,shadow}.
\]

In parallelo viene ricostruita la deformazione differenziale della wing tip.
La sottrazione elimina, al primo ordine, il moto rigido e la risposta forzata
della manovra. Lo stesso filtro band-pass, lo stesso inviluppo di Hilbert e la
stessa regressione logaritmica sono usati in ogni campagna.

## Manovra

La sequenza è:

1. volo iniziale controllato;
2. inizio manovra a 8 s;
3. rampa di picchiata, limitata a 2°/s e mai più breve di 2 s;
4. hold in picchiata per 1 s;
5. richiamata lineare di 4.05 s;
6. SAS-off due secondi dopo l'inizio della richiamata, per gli ultimi 2.05 s;
7. hold di tutte le superfici durante SAS-off, con il solo rap noto nella
   excited;
8. fine della run di produzione due passi dopo la chiusura della finestra
   SAS-off.

Il generatore sorgente definisce anche ritorno al riferimento e recupero, già
osservati nella campagna precursore. Non vengono integrati nelle nuove run:
avvengono dopo tutti i campioni usati per stimare il polo e non possono quindi
modificare il risultato. Il tratto iniziale fino al rilascio è invece conservato
integralmente, perché serve a stabilizzare gli stati di Wagner, i PID, i filtri
e gli attuatori.

La classe nominale non impone `n` in feedback. Definisce il pitch-rate

\[
q_c=\frac{(n_{nom}-1)g}{V}
\]

e l'ampiezza centrata

\[
A=\frac{1}{2}q_c(2.00+2.05).
\]

Il fattore di carico è sempre ricostruito dalla shadow tramite accelerazione
specifica nel frame body. Sono salvati media, deviazione standard, estremi e
pendenza di (n), oltre a pitch-rate e TAS nella finestra SAS-off.

## Ottimizzazione del costo

La campagna iniziale richiedeva 78 traiettorie e 235443 passi. Quella pronta ne
richiede 28 e 49906 passi: **-64% di traiettorie e -79% di integrazione**. Le
riduzioni sono fisicamente motivate:

- i risultati precursori a 65, 66 e 67 m/s collocano tutti i crossing vicino a
  66.75 m/s; bastano quindi tre punti 66.25/66.75/67.25 per bracket, curvatura
  locale e interpolazione;
- tre livelli 1.0/1.3/1.6 permettono di vedere una tendenza e la sua eventuale
  non linearità. `n=1.8` è escluso perché nei risultati esistenti la deviazione
  standard della TAS è 0.56–0.74 m/s, oltre il limite 0.50 m/s: quelle coppie
  non sarebbero utilizzabili come poli locali quasi stazionari;
- `timestep` genera solo il caso fine: il riferimento a 0.01 s è lo stesso caso
  già prodotto da `primary`;
- `stiffness_screen` genera solo -1% e +1%: il punto centrale è riusato da
  `primary`. Tre punti centrati sono sufficienti per la derivata locale;
- ogni run termina dopo la finestra identificata; il recupero non osservato è
  eliminato;
- il NetCDF campiona a 0.01 s solo dagli ultimi 4 s prima del rilascio e salva
  esclusivamente i sette nodi e il joint modale letti dall'analisi. Il passo di
  integrazione a 0.005 s resta realmente fine, anche se l'output è a 0.01 s.

Non è stato usato il restart binario per dividere shadow/excited: il manuale
MBDyn locale lo dichiara sperimentale e sostanzialmente abbandonato. Con PID,
filtri discreti e sample-and-hold sarebbe un rischio di riproducibilità maggiore
del risparmio residuo.

## Matrice delle simulazioni

### `primary`: risultato fisico principale

- TAS: 66.25, 66.75 e 67.25 m/s;
- classi nominali: 1.0, 1.3 e 1.6;
- passo: 0.01 s;
- rigidezza: FEM lineare originale;
- 9 condizioni × shadow/excited = **18 traiettorie**.

La griglia è centrata sull'onset paired precedente, circa 66.75 m/s, e conserva
0.50 m/s su entrambi i lati. Il caso 1 g usa lo stesso metodo delle manovre: non
viene confrontato con poli SSI ottenuti da un'altra procedura.

### `timestep`: convergenza numerica

- TAS: 66.75 m/s;
- classe: 1.6;
- passo generato: 0.005 s;
- riferimento riusato da `primary`: 0.01 s;
- 1 condizione × shadow/excited = **2 nuove traiettorie**.

Cambiare soltanto `TIME_STEP` sarebbe scorretto, perché i filtri discreti e gli
attuatori cambierebbero dinamica fisica. Il generatore ridiscretizza con Tustin
prewarped:

- Butterworth laterale del secondo ordine: 0.6 Hz;
- filtro quota: 0.15 Hz;
- filtri velocità verticale e pitch-rate: 1 Hz;
- attuatore del primo ordine: \(\tau=0.01\) s.

Per esempio, a 0.005 s l'attuatore diventa
`A1=0.6, B0=B1=0.2`, preservando lo stesso sistema continuo.

### `stiffness_screen`: gate prima del prestress

- TAS: 66.75 m/s;
- classi: 1.0 e 1.6;
- scale generate del modo 7: 0.99 e 1.01;
- scala 1.00 riusata da `primary`;
- passo: 0.01 s;
- 4 condizioni × shadow/excited = **8 nuove traiettorie**.

Poiché la massa modale è unitaria e `K_77=408.5998`, una scala di frequenza
`s_f` richiede

\[
\Delta K_{77}=(s_f^2-1)K_{77}.
\]

Il caso MBDyn aggiunge la forza generalizzata

\[
Q_7^{screen}=-\Delta K_{77}q_7,
\]

con rampa 0–1 nei primi 5 s, così il controllo riassesta l'aeromobile prima
della manovra. Viene perturbato solo il modo dominante e non cambiano forme
modali o accoppiamenti. È quindi una derivata parametrica locale, **non un
prestress fisico**.

L'analisi interpola \(d\sigma/ds_f\) e valuta lo span previsto di \(\sigma\) per
uno spostamento ±1% della frequenza. Se supera 0.05 1/s, raccomanda il ROM
fisico; in caso contrario il prestress non è risolvibile con la finestra e il
metodo attuali e non serve per la conclusione principale.

## Criteri di qualità

Una coppia può contribuire all'onset soltanto se:

- errore tra `n` medio e classe nominale ≤ 0.20;
- deviazione standard di `n` ≤ 0.18;
- pendenza assoluta di `n` ≤ 0.30 1/s;
- errore fra TAS media e TAS nominale ≤ 1.00 m/s;
- deviazione standard della TAS ≤ 0.50 m/s;
- errore del pitch-rate medio ≤ max(0.75°/s, 30% del comando);
- deviazione standard del pitch-rate ≤ 1.50°/s;
- differenza tra \(\sigma_{q7}\) e \(\sigma_{tip}\) ≤ 0.25 1/s;
- contaminazione differenziale delle superfici nella banda BFF ≤ 0.025°.

La TAS è calcolata come norma della velocità dell'aria relativa al velivolo,
`[V_INF,0,0] - XP_base`; `XP_base` da solo è una velocità perturbativa, non la
TAS. Le soglie di stazionarietà sono applicate nella finestra di identificazione.
Quota e recupero non entrano nel polo open-loop.

Per la convergenza temporale sono richiesti:

- \(|\sigma_{0.005}-\sigma_{0.01}|\le0.05\) 1/s;
- \(|n_{0.005}-n_{0.01}|\le0.03\).

La zona \(|\sigma|\le0.05\) 1/s è classificata `near_onset`. L'onset viene
interpolato linearmente soltanto fra due velocità valide con cambio di segno.

## Come eseguire

### 1. Verifica prima delle run

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 audit_modal_joint.py
python3 verify_setup.py
```

`verify_setup.py` controlla numero di input, copertura paired, metadata,
`TIME_STEP`, conteggio delle forze e hold delle superfici. Il report è
`audit/preflight.json`.

### 2. Campagna principale

```bash
python3 run_sweep.py --campaign primary --execute --jobs 2 --analyse
```

### 3. Convergenza temporale

```bash
python3 run_sweep.py --campaign timestep --execute --jobs 2 --analyse
```

L'analisi aggiunge automaticamente `BFF_PULLUP_V2/primary` come riferimento a
0.01 s.

### 4. Screening di rigidezza

```bash
python3 run_sweep.py --campaign stiffness_screen --execute --jobs 2 --analyse
```

Anche qui l'analisi riusa automaticamente i punti centrali di `primary`.
L'ordine sopra è quindi obbligatorio. Se analizzati manualmente, i comandi
equivalenti sono:

```bash
python3 analyse_sweep.py /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/timestep \
  --reference-directory /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/primary
python3 analyse_sweep.py /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/stiffness_screen \
  --reference-directory /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/primary
```

`--jobs 2` è la scelta consigliata sulla macchina attuale; `--jobs 1` resta la
più conservativa. In base ai tempi reali precedenti, la nuova campagna richiede
circa **4.8–5.5 ore seriali** o **2.7–3.3 ore con due processi**, contro circa
22.5 ore seriali della matrice iniziale. La stima non accredita alcun guadagno
di CPU dalla riduzione dell'output, quindi è prudente. Lo spazio atteso è
dell'ordine di decine di MB e va verificato sul primo caso, invece degli 8–10 GB
stimati con l'output aerodinamico completo.

Senza `--execute`, `run_sweep.py` rigenera input e manifest ma non lancia MBDyn.
Un NetCDF viene riutilizzato solo se raggiunge davvero la fine di SAS-off;
`--overwrite` sostituisce soltanto gli output del caso con lo stesso stem.

Per cambiare la matrice, modificare `campaign_config.json` e rigenerare gli
input. I manifest JSON contengono gli SHA-256 delle sorgenti da cui gli input
sono stati costruiti.

## Output dell'analisi

Ogni `C:\Users\Utente\Desktop\BFF_PULLUP_V2\<campagna>\analysis` contiene:

- `paired_results.csv`: stato di volo, \(\sigma\), controlli qualità e verdict;
- `onset.csv`: bracket e interpolazione per classe, passo e scala;
- `timestep_convergence.csv`: confronto 0.01/0.005 s, quando disponibile;
- `stiffness_sensitivity.csv`: regressione \(\sigma\)-frequenza, quando disponibile;
- `summary.json`: riepilogo completo e gate prestress;
- `summary.png`: mappe sintetiche dei risultati validi.

L'analisi rifiuta run incomplete e coppie mancanti. Non sostituisce una shadow
di una condizione con quella di un'altra; riusa tra campagne soltanto una coppia
fisicamente identica, dichiarata come riferimento.

## Risultati ottenuti il 2 settembre 2026

Le 28 traiettorie sono terminate senza fallimenti: 18/18 `primary`, 2/2
`timestep` e 8/8 `stiffness_screen`. Tutti i NetCDF raggiungono la fine della
finestra SAS-off, non ci sono coppie mancanti e tutti i punti superano sia i
criteri di stazionarietà sia quelli di identificazione.

### Onset con rigidezza lineare originale

| Classe nominale | `n` medio ottenuto | Bracket [m/s] | Onset interpolato [m/s] |
|---:|---:|---:|---:|
| 1.0 | 0.9983 | 66.75–67.25 | 66.9022 |
| 1.3 | 1.3099 | 66.75–67.25 | 66.8341 |
| 1.6 | 1.6395 | 66.75–67.25 | 66.8098 |

La tendenza è una riduzione dell'onset di 0.068 m/s a `n=1.3` e 0.092 m/s a
`n=1.6` rispetto a 1 g. È piccola rispetto alla spaziatura di 0.50 m/s della
griglia: va riportata come tendenza interpolata del ROM lineare, non come
spostamento sperimentale risolto con incertezza inferiore a 0.1 m/s.

### Convergenza temporale

Nel caso 66.75 m/s, `n=1.6`:

- `sigma(dt=0.01) = -0.03617 1/s`;
- `sigma(dt=0.005) = -0.02008 1/s`;
- differenza `+0.01610 1/s`, inferiore alla soglia 0.05 1/s;
- differenza nel fattore di carico `-0.00115`, inferiore alla soglia 0.03.

La verifica è quindi `converged=True`: il passo 0.01 s è adeguato alla
risoluzione dichiarata dello studio.

### Screening di rigidezza

| Classe | dσ per +1% di frequenza [1/s] | Span σ previsto fra -1% e +1% [1/s] | RMS fit [1/s] |
|---:|---:|---:|---:|
| 1.0 | +0.11192 | 0.22384 | 0.04130 |
| 1.6 | +0.12428 | 0.24856 | 0.04006 |

Entrambi gli span sono molto superiori alla risoluzione 0.05 1/s. Il gate
finale è pertanto `physical_stress_stiffening_rom_recommended`. Nel modello
parametrico, a 66.75 m/s, -1% resta stabile mentre +1% diventa instabile; questo
segno non va però attribuito direttamente alla richiamata. Lo screening misura
la derivata rispetto a una variazione imposta di `K77`, mentre il modello
attuale non determina né segno né ampiezza della variazione fisica prodotta dal
campo di sforzo.

Conclusione: il risultato lineare sull'onset è numericamente convergente, ma
per affermare se il prestress della richiamata anticipa o ritarda il flutter
serve ora il ROM fisico descritto nella sezione successiva e in
`nastran/README.md`.

## Verifica fisica frozen-time preparata ed eseguita lato MBDyn

Per non spendere ore in una nuova matrice sono state aggiunte soltanto due
shadow a 66.75 m/s, `n_nom=1.0` e `n_nom=1.6`, senza rap BFF. Conservano
manovra, controllo, passo 0.01 s e modello delle run primarie, ma salvano i 58
elementi aerodinamici, tre punti di Gauss per elemento, il nodo base e i 45
nodi di interfaccia. L'output resta limitato agli ultimi secondi e ha passo
effettivo 0.02 s (50 Hz).

`extract_prestress_loads.py` usa la parte interna della finestra SAS-off,
scartando 0.25 s a ciascun estremo:

| Stato | Finestra [s] | Campioni | n corpo medio | std(n) | Fz aero [lbf] | std(Fz) [lbf] |
|---|---:|---:|---:|---:|---:|---:|
| 1 g | 13.25–14.79 | 78 | 0.9940 | 0.0020 | 417.278 | 0.878 |
| 1.6 g | 16.38–17.90 | 77 | 1.6329 | 0.0356 | 694.622 | 21.646 |

Per ogni campione le forze e i momenti sono trasformati dal riferimento globale
al body frame con la matrice del nodo modale. Le risultanti per unità di apertura
sono integrate con i pesi di Gauss `5/9, 8/9, 5/9` e trasportate al relativo
reference grid RBE3. I corpi delle superfici mobili 880xxx/881xxx sono riportati
al corrispondente nodo strutturale 990xxx/991xxx.

Il campo inerziale Nastran contiene

\[
F_{inerziale}=M\,(g-a_{CG}),
\]

tramite `GRAV`, più `RFORCE` attorno al grid fusoliera 10062 per velocità
angolare e accelerazione angolare opposta (carico di d'Alembert). Il grid 10062
è a meno di 0.8 in dal CG del modello. Lo stesso grid è vincolato nel subcase
statico e nello screening modale, così baseline scarica e stati precaricati
hanno identiche condizioni al contorno. Questo evita di attribuire al prestress
lo scarto prodotto dal passaggio free-free/supportato; i modi vengono poi
identificati rispetto alla base MBDyn free-free mediante MAC.

Prima del bilanciamento il residuo verticale vale 0.364 lbf a 1 g e 9.718 lbf a
1.6 g. Il residuo totale in forza è 3.96% e 6.54% della somma delle norme ed è
quasi tutto assiale: 33.22 e 90.31 lbf. È la reazione di trim/thrust implicita
del `total pin joint` MBDyn, che blocca le traslazioni X/Y. La corrispondente
forza e coppia sono applicate al grid già vincolato: chiudono le risultanti senza
alterare il campo di tensione della parte libera. Il residuo in momento è molto
piccolo, rispettivamente 0.285 e 1.003 lbf-in in norma.

Il contributo DLM non ha una distribuzione nodale unica perché MBDyn lo applica
come sola forza generalizzata del modo 7. Per rendere la decisione robusta sono
forniti due deck per ogni stato:

- `nodlm`: carichi fisici c81/Theodorsen e inerziali;
- `with_dlm`: aggiunge la realizzazione canonica
  `M_lumped*phi_7*Q_7/(phi_7^T*M_lumped*phi_7)` su tutti gli 8527 grid.

La proiezione della seconda distribuzione sul modo 7 è esattamente 1. Il suo
valore medio è `Q7=+0.7951` a 1 g e `Q7=-60.1338` a 1.6 g. I due risultati
costituiscono un bracket d'incertezza; la versione `with_dlm` non va descritta
come una ricostruzione dei carichi di pannello DLM.

I deck sono in:

```text
C:\Users\Utente\Desktop\BFF_PULLUP_V2\load_recovery\prestress_loads\n1p0
C:\Users\Utente\Desktop\BFF_PULLUP_V2\load_recovery\prestress_loads\n1p6
```

Da WSL, con `nastran` nel `PATH`:

```bash
cd /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/load_recovery/prestress_loads/n1p0
./run_nastran.sh
cd ../n1p6
./run_nastran.sh
```

Se il comando ha un altro nome:

```bash
NASTRAN_CMD=/percorso/al/comando/nastran ./run_nastran.sh
```

Su Zeno i deck sincronizzati sono in `~/ZENO/prestress_loads/n1p0` e
`~/ZENO/prestress_loads/n1p6`. Il wrapper locale richiede l'ordine diverso
`nast old=no deck.bdf`; dalla cartella padre usare quindi:

```bash
cd ~/ZENO/prestress_loads
./run_all_zeno.sh
```

Il generatore non scrive card `FORCE` o `MOMENT` con vettore direzione nullo;
il controllo è necessario per evitare `USER FATAL 9994 (BULFUN)` nella
distribuzione modale DLM.

La decomposizione prestressata usa `PARAM,TESTNEG,3`, come definito dalla Quick
Reference MSC, per estrarre e diagnosticare anche eventuali radici negative.
Il postprocess le riporta esplicitamente e blocca il ROM se ne è presente una:
il parametro permette di misurare l'instabilità, non di accettarla.

La cartella `n1p0` include inoltre `supported_unloaded.bdf`, una singola SOL 103
scarica con lo stesso supporto `SPC=620`. Questo riferimento elimina dal
confronto lo scarto dovuto alle diverse condizioni al contorno tra il FEM
free-free originale e lo screening prestressato supportato al CG.

Poi:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 analyse_prestressed_modes.py
```

Il postprocess abbina i modi FEM 7–12 con MAC sui 45 nodi di interfaccia e
confronta la variazione 1 g→1.6 g con la soglia risolvibile dello 0.4%. Il ROM
stress-stiffened va costruito soltanto se entrambi i bracket hanno MAC almeno
0.90, superano la soglia e concordano sul segno. La procedura completa e i
controlli da fare nei `.f06` sono in `nastran/README.md`.

Il percorso dinamico finale preferito resta il `RECORD GROUP 19`: i SOL 103
preparati sono una verifica frozen-time e non aggiornano `Kgeo(t)` durante tutta
la richiamata.

## Risultati dello screening prestress e rettifica del gate

Le cinque SOL 103 eseguite su Zeno (baseline supportata e quattro combinazioni
1 g/1.6 g, senza/con DLM) terminano con `END OF JOB`. Con
`PARAM,TESTNEG,3` sono state estratte e controllate anche le possibili radici
negative, ma il postprocess precedente scartava quelle sotto 0.1 Hz in modulo
di frequenza equivalente. Questa soglia non dimostra che siano numeriche.
A 1.6 g sono presenti λ=-0.03927591 s⁻² senza DLM e λ=-0.07955200 s⁻² con
DLM: la stabilità fisica non è accertata. Occorre diagnosticare le forme e
il vincolo, non interpretare `END OF JOB` come certificazione di stabilità.

Il modo FEM 7 è associato al modo supportato 4 con `MAC=0.9197–0.91975`. A
parità di `SPC=620` la frequenza scarica è 2.768387 Hz. I risultati sono:

| Stato | senza DLM [Hz] | shift | con DLM [Hz] | shift |
|---|---:|---:|---:|---:|
| 1 g | 2.763266 | -0.1850% | 2.763275 | -0.1847% |
| 1.6 g | 2.752805 | -0.5629% | 2.751967 | -0.5931% |

Lo shift incrementale 1 g→1.6 g è -0.3786% senza DLM e -0.4092% con DLM.
Le due realizzazioni concordano sul softening nel problema supportato. Tuttavia
il ramo con DLM usa un precarico ricostruito con coefficienti errati, e la
distribuzione nodale convenzionale non delimita rigorosamente il vero carico
DLM. Il precedente codice confrontava lo shift dalla baseline scarica con
0.4%, non necessariamente l'incremento 1 g→1.6 g. La decisione automatica
`physical_prestress_resolved_robust_to_dlm_distribution_build_rom` va quindi
considerata non validata, non una conclusione fisica acquisita.

I modi 8, 10, 11 e 12 cambiano meno dello 0.18%; il matching del modo 9 ha MAC
circa 0.814 e non viene usato per la decisione. La sensibilità MBDyn ±1% indica
che lo shift del modo 7 può produrre una variazione del tasso di crescita di
ordine 0.07 1/s, maggiore della differenza di convergenza temporale osservata
(0.0161 1/s). È una motivazione per un confronto mirato, non una misura già
validata dell'effetto fisico del prestress sull'onset free-free.

Prima di validare il ROM vanno risolte le anomalie riportate. I modi SOL 103 supportati non
devono essere inseriti direttamente nel modal joint free-free. La realizzazione
corretta richiede le matrici geometriche ridotte nel `RECORD GROUP 19`, oppure
un'esportazione FEM equivalente che conservi base e condizioni free-free.

## File creati e controlli effettuati

- `campaign_config.json`: unica matrice di progetto e soglie;
- `campaign.py`: rendering, comando di manovra, ridiscretizzazione e forza di
  screening;
- `run_sweep.py`: preparazione e lancio esplicito con manifest incrementale;
- `analyse_sweep.py`: identificazione paired, qualità, onset, convergenza e
  gate prestress;
- `audit_modal_joint.py`: audit FEM/manuale/sorgente MBDyn;
- `verify_setup.py`: preflight statico riproducibile;
- `run_load_recovery.py`: due sole shadow con output di recupero carichi;
- `extract_prestress_loads.py`: integrazione aero, inerzia, equilibrio e deck
  SOL 103 con/senza distribuzione DLM;
- `analyse_prestressed_modes.py`: tracking MAC dei modi 7–12 e gate fisico;
- `C:\Users\Utente\Desktop\BFF_PULLUP_V2\*\cases`: 28 input MBDyn
  ottimizzati e tutti i risultati prodotti dalle run;
- `runs_superseded_full_matrix`: copia conservativa della matrice estesa, da
  non eseguire;
- `audit/modal_joint_audit.json` e `audit/preflight.json`: evidenze dell'audit;
- `C:\Users\Utente\Desktop\BFF_PULLUP_V2\load_recovery`: NetCDF di recupero,
  riepiloghi, carichi nodali e quattro deck Nastran pronti;
- `results/prestress`: CSV e JSON conclusivi del gate fisico;
- `nastran/README.md`: procedura, risultati e criteri Nastran.

Sono stati eseguiti `py_compile`, audit, preflight, parsing MBDyn e tutte le 30
integrazioni MBDyn: 18 primary, 2 timestep, 8 stiffness screen e 2 load
recovery. Tutti i controlli hanno avuto esito positivo. È stata inoltre
verificata la presenza in entrambi i NetCDF di `X`, `F` e `M` per tutti i 58
elementi e l'assenza di carte Bulk dati oltre 80 colonne. Su Zeno sono state
eseguite con successo la baseline supportata e le quattro run prestressate.

## Limiti interpretativi

Il modello aerodinamico resta la combinazione Theodorsen/c81 più correzione
DLM calibrata sul punto 7. Lo screening varia una sola diagonale modale e non
predice il campo di sforzo reale. La classe `n=1.8` non è inclusa perché la
campagna precursore la mostra non quasi-stazionaria con la soglia adottata. La
finestra SAS-off è breve: \(\sigma\) è un
valore locale medio se `n`, TAS o pitch-rate variano. Per questo i risultati
vanno sempre letti insieme agli indicatori di stazionarietà e alla coerenza
modo 7–wing tip.

La SOL 103 usa un singolo campo medio frozen-time: non rappresenta la variazione
continua del prestress. Il thrust assiale non è distribuito sui propulsori,
perché anche il modello MBDyn lo rappresenta come reazione del pin al CG. Il
confronto senza/con DLM è una sensibilità a una realizzazione generalizzata,
non un bracket certificato, e non sostituisce il recupero di carichi di pannello.
Va inoltre corretto il recupero dei coefficienti dai deck. Questi limiti devono restare
espliciti nell'interpretazione del gate fisico.
