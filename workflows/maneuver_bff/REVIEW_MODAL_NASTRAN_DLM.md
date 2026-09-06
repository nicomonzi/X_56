# Modal joint, Nastran, DLM e confronto lungo la baseline

Revisione del 5 settembre 2026. Analisi del codice e dei risultati esistenti;
nessuna nuova integrazione MBDyn, nessuna run Nastran e nessuna modifica ai deck.
Questa revisione corregge le conclusioni del precedente gate sul prestress.

## Conclusione operativa

Il modello corrente rappresenta gli effetti inerziali e aerodinamici della
manovra, ma non la variazione di rigidezza elastica causata dagli sforzi.
MBDyn può rappresentarla. Mancano però le matrici geometriche nella stessa
base free-free impiegata dal modal joint. Le frequenze SOL 103 attualmente
disponibili non bastano a ricostruirle.

Consiglio un primo confronto della dinamica perturbata lungo la baseline
esistente, con una matrice geometrica 6×6 dipendente dal tempo, ricostruita
offline da pochi stati/carichi indipendenti. Non consiglio un altro sweep
completo, né di sostituire direttamente i modi supportati Nastran a quelli
free-free. Prima vanno risolte due anomalie concrete: radici negative filtrate
dal postprocess e coefficienti DLM sbagliati nel recupero del precarico.

## 1. Cosa fa realmente il modal joint

Il modello usa il modal joint 5, il nodo di riferimento 990000 e 45 nodi di
interfaccia. Il FEM contiene 8527 nodi e 60 modi; l'integrazione seleziona i
sei modi elastici 7–12, mass-normalizzati. Le frequenze sono 3.217134,
5.302727, 8.705123, 11.163976, 12.257061 e 12.758920 Hz.
Lo smorzamento strutturale prescritto è 1%.

La cinematica floating-frame è

\[
x_P=x_0+R_0(f_P+\Phi_Pq).
\]

Sono ammesse grandi traslazioni e rotazioni del riferimento; la deformazione
elastica resta una combinazione lineare di forme modali e deve rimanere piccola.
La dinamica del sistema complessivo non è per questo lineare: cambiano assetto,
forze aerodinamiche, velocità locali, accoppiamenti inerziali e controlli.
Non è però presente una legge di rigidezza strutturale dipendente dal precarico.

Il file contiene i gruppi 1–11, non il gruppo 19. Il gruppo 11 fornisce
l'inerzia nodale con cui vengono ricostruiti gli invarianti inerziali: non è
una matrice geometrica. Non bisogna confondere stress-stiffening con termini
centrifughi, Coriolis o spin-softening già presenti nella formulazione inerziale.

### Il percorso nativo: RECORD GROUP 19

Il manuale tecnico e `modalad.cc` implementano una correzione quasi-statica:

\[
K_{geo}(t)=\sum_j c_j(t)G_j,
\qquad \delta W\simeq\delta q^T[K_0+K_{geo}(t)]q.
\]

Ogni matrice costante `G_j` deriva dal campo di sforzo prodotto da un carico
unitario. Il relativo coefficiente viene aggiornato durante l'integrazione.

| Tag nel gruppo 19 | Coefficiente dinamico, espresso nel riferimento corpo |
|---|---|
| 1–6 | ωx², ωy², ωz², ωxωy, ωyωz, ωxωz |
| 7–9 | accelerazione angolare |
| 10–12 | accelerazione traslazionale meno gravità |
| 13–15, 16–18 | forza e momento del primo nodo di interfaccia |
| 19–21, 22–24, … | forza e momento delle interfacce successive |

Nel codice locale i tag delle forze fanno riferimento all'ordine dei nodi
di interfaccia del joint, non all'etichetta numerica del GRID.
Con 45 interfacce ci sono al massimo 282 canali; si può fornire un sottoinsieme.
Omettere un canale significa assumerne nullo il contributo: la scelta va validata.
Le matrici del file devono essere coerenti con `NMODES=60`; il lettore estrae
poi la sottobase attiva. Non si possono appendere matrici 6×6 a un header da
60 modi senza una trasformazione coerente del formato.

Nella revisione locale il percorso stress-stiffening è quello `ModalAd`,
selezionato abilitando automatic differentiation; il modello corrente non
lo abilita. La `femgen` locale non produce il gruppo 19. Occorrono quindi
esportazione/riduzione FEM e verifica del nuovo input, non una singola opzione.

**Limite importante:** le forze del gruppo 19 sono ricavate dalle reazioni
del joint alle interfacce. Una forza applicata direttamente con `force: ...,
modal`, come la correzione DLM o il damper modale, non genera automaticamente
quel canale di recupero dello sforzo. Per includerne il precarico serve una
realizzazione fisica esplicita, oppure un contributo geometrico separato
associato alla forza generalizzata. Non basta attivare il gruppo 19.

Fonti locali: `src/mbdyn/manual/tecman/tecman-modal.tex`, sezione
“Quasi static corrections for stress stiffening”; `manual/input/modal-fem-format.tex`;
`mbdyn/struct/modal.cc`, `modalad.cc`, `modalforce.cc` sotto `/home/nicomonzi/src/mbdyn`.

## 2. Aerodinamica attuale e significato della correzione DLM

### Parte distribuita sulle sezioni

`BFF_open_loop/INCLUDE/aerobody.mbd` contiene 58 elementi aerodinamici con
tre punti di Gauss ciascuno: 174 sezioni. Ogni sezione usa la velocità locale
del nodo deformato/ruotato, la velocità angolare e il vento. Il modello
Theodorsen realizza il ritardo circolatorio con due stati per sezione
(348 stati aerodinamici complessivi), più contributi non circolatori.

I parametri di default locali sono A1=0.165, A2=0.335, b1=0.0455, b2=0.3;
a velocità costante corrispondono all'approssimazione di Wagner
`1 - A1 exp(-b1 τ) - A2 exp(-b2 τ)`, con `τ=2Ut/c`.
La collocazione della condizione aerodinamica è spostata di mezza corda dal
punto di forza: forza al quarto di corda e condizione ai tre quarti.
Questo non è un calcolo DLM tridimensionale istantaneo e non risolve una
scia tridimensionale/pannellatura globale nella manovra.

La tabella `x56_effective.c81` deriva da NACA 0012, ma è una polare efficace
calibrata: CL è scalato di 0.7732; la pendenza aggiunta a Cm è -0.005930 per
grado con raccordo alle alte incidenze; CD rimane quello tabulato di partenza.
Il report precedente documenta l'accordo delle pendenze globali con Nastran:
CZα=0.108331/deg contro 0.108327/deg, Cmα=-0.007550/deg contro -0.007563/deg.
È una calibrazione dell'aeroplano, non una nuova misura sperimentale di profilo.

Il recupero dei carichi ai punti di Gauss è strutturalmente coerente con il
codice `aeroelem.cc`: `F_gp` e `M_gp` sono per unità di apertura, prima della
moltiplicazione per peso di Gauss e semispan. `X_gp` è il riferimento della
sezione; il momento va trasportato al nodo con `(X_gp-X_node)×F_gp`.
Non va aggiunto una seconda volta lo stesso braccio aerodinamico di corda.

### Parte calibrata sui risultati SOL 145

Nel deck la correzione è applicata direttamente al solo modo 7:

\[
Q_{DLM,7}=r(t)[K_7(V_{nom})(q_7-q_{eq}(V_{nom}))+C_7(V_{nom})\dot q_7],
\qquad r(t)=\min(1,t/5).
\]

`run_case.py` interpola una tabella di coefficienti tra 50 e 70 m/s.
Questi sono fissati al rendering sul valore nominale della TAS: non si
aggiornano con la TAS istantanea durante la manovra. Il termine resta attivo
anche quando il SAS viene spento. I coefficienti positivi sul membro delle
forze riducono rigidezza e smorzamento efficaci del modo, non li aumentano.

Questo codice implementa una **correzione modale di rango uno calibrata sui
poli SOL 145**, non una RFA/Roger ricavata esplicitamente dalle matrici DLM
multimodali `Qhh(k,M)`. La denominazione “DLM ROM” nei file va letta in questo
senso limitato. Non è giustificato dedurre da questi due coefficienti un campo
univoco di pressioni/sforzi.

L'accordo non è perfetto: il report di baseline riporta a 65 m/s frequenza
MBDyn 1.9980 Hz contro 2.0552 Hz e σ=-0.4435 contro -0.2057 s⁻¹; a 70 m/s
1.6852 contro 2.0644 Hz. Il riferimento SOL 145 raffinato ha onset 65.7637 m/s
TAS, mentre lo screening di manovra ha un onset locale a 1 g di circa 66.9022.
Il confronto di differenze sulla stessa baseline è quindi più difendibile
dell'affermazione di un onset assoluto accurato al centesimo di m/s.
Per convertire i risultati Nastran va mantenuta la distinzione KEAS/TAS
documentata dal parser; non basta confrontare direttamente la colonna VELOCITY.

Durante lo studio di rigidezza, lascerei invariati questi coefficienti:
ricalibrarli dopo ogni modifica strutturale confonderebbe effetto strutturale
ed effetto della nuova taratura. Se in futuro serve un vero modello DLM
multimodale, occorrono matrici generalizzate comprensive dei moti rigidi e
dei comandi; aggiungerle integralmente alla strip theory conterebbe due volte
parte dell'aerodinamica. Bisogna sostituire la parte lineare o identificare
una correzione residua coerente.

## 3. Modello Nastran e cosa dicono davvero le run esistenti

È un modello misto, non una singola trave: il riepilogo F06 conta 8527 GRID,
1682 CBAR, 222 CBEAM, 181 CBUSH, 316 CHEXA, 158 CPENTA, migliaia di shell,
127 PCOMP, 639 CONM2, 61 RBE3 e 7 RBE2. Le unità del flusso sono in-lbf-s
con `PARAM,WTMASS,0.002591`; coordinate e trasformazioni dei GRID/MPC devono
restare coerenti durante qualsiasi proiezione.

Le run in `/home/nicomonzi/ZENO/prestress_loads` fanno una statica di
precarico seguita da SOL 103 con `STATSUB=1`, usando `SPC=620` su tutti i sei
DOF del GRID 10062 anche nella modale. Il GRID è vicino, ma non coincidente,
al CG. Sono quindi modi **supportati**, diversi dal modello free-free del
modal joint. Il riferimento scarico con gli stessi SPC è utile per separare
l'effetto di carico da quello del vincolo, ma non elimina la necessità di
verificare l'influenza del vincolo sul campo di sforzo e sulla sottobase.

La SOL 103 prestressata usa rigidezza lineare più differenziale; quest'ultima
può comprendere anche contributi follower. Nell'esportazione destinata a MBDyn
va separata la parte di stress geometrico dai termini già rappresentati
dall'aerodinamica e dall'inerzia, per evitare duplicazioni.
[Fonte: MSC Nastran Dynamic Analysis User's Guide, “Prestiffened Normal Mode Analysis”](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2023.2_Dynamic_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2023.2_Dynamic_Analysis_User_Guide.pdf).

Il recupero usa le shadow a V=66.75 m/s nelle finestre interne SAS-off:
13.25–14.79 s a 1 g, 16.38–17.90 s a 1.6 g. Il carico specifico verticale
medio misurato è rispettivamente 0.99398 g e 1.63290 g; nel secondo caso la
deviazione standard è 0.03555 g. Quindi `nnom=1.6` non è esattamente 1.6 g
costante. Vanno conservati gravità reale, accelerazioni e velocità angolari,
non sostituiti da una gravità fittizia moltiplicata per n.

### Anomalia A: radici negative effettivamente presenti

Il postprocess precedente esclude gli autovalori negativi con modulo minore
di `(2π·0.1)^2`. Questa soglia non dimostra che siano numerici.

| Caso | Prima radice λ [s⁻²] | Colonna frequenza F06 [Hz] |
|---|---:|---:|
| Supportato scarico | +0.3772218 | 0.09775039 |
| 1 g senza DLM | +0.3114289 | 0.08881771 |
| 1 g con DLM | +0.3116980 | 0.08885607 |
| 1.6 g senza DLM | **-0.03927591** | 0.03154157 |
| 1.6 g con DLM | **-0.07955200** | 0.04488960 |

La frequenza positiva stampata non rende positivo l'autovalore. Le due run
a 1.6 g hanno una direzione negativa della rigidezza combinata nel problema
calcolato. Resta da stabilire se sia una reale instabilità o un problema di
vincolo/connettività/formulazione: servono la forma completa della radice,
i DOF presso il supporto, energie e controlli di vincolo. `END OF JOB` e
`TESTNEG=3` non certificano la stabilità. Non si devono azzerare le radici o
ignorare il fatal per trasformare il problema in un caso fisicamente accettato.

### Anomalia B: carico DLM recuperato con coefficienti diversi dal deck

`extract_prestress_loads.py:241` usa valori costanti non corrispondenti ai
deck realmente eseguiti. Il ricalcolo sui NetCDF originali dà:

| Quantità | Recupero precedente | Deck realmente eseguito |
|---|---:|---:|
| K7 | 146.378920 | 147.100000 |
| C7 | 3.678940 | 10.381000 |
| q7 equilibrio | -0.4304727364 | -0.41888588052 |
| Q7 media, 1 g | +0.795083 | **-0.891132** |
| Q7 media, 1.6 g | -60.133822 | **-62.610336** |

Le Q7 sono forze generalizzate nella normalizzazione modale del modello,
non singole forze nodali in lbf. A 1.6 g l'errore di modulo è circa 4%; a 1 g,
vicino alla compensazione del trim, cambia anche il segno. Questa incoerenza
non cambia le traiettorie MBDyn già integrate: interessa la successiva
ricostruzione del precarico dei deck `with_dlm`.

Inoltre, la realizzazione `M_lumped φ7 Q7/(φ7ᵀ M_lumped φ7)` è una scelta
convenzionale. Infinite distribuzioni hanno la stessa proiezione modale e
campi di sforzo diversi. Vanno verificati anche la proiezione sui modi 8–12
e il risultante rigido dopo il bilanciamento. Il confronto senza/con questa
distribuzione è una sensibilità a due ipotesi, **non un intervallo certificato**
che contenga il vero effetto delle pressioni DLM.

### Cosa rimane utilizzabile

Il ramo associato al modo free-free 7, individuato come modo supportato 4,
scende da 2.763266 a 2.752805 Hz tra 1 e 1.6 g senza carico DLM ricostruito:
-0.3786%. La variante con il recupero precedente scende del -0.4092%.
È un'indicazione di **softening nel problema supportato**, non una prova che
la vera struttura free-free irrigidisca o ammorbidisca della stessa quantità.
Stress-stiffening è il nome del termine: il suo effetto può avere entrambi i segni.

Il MAC precedente mescola traslazioni in inches e rotazioni in radianti e usa
solo 45 interfacce. Va sostituito/affiancato con MAC mass-weighted oppure
traslazionale e rotazionale separati con scaling esplicito, confrontando
prima le forme supportate fra loro. Il modo 9 ha inoltre MAC circa 0.814.
Non è lecito trasformare il gate automatico precedente in una validazione ROM.

## 4. Trasferimento Nastran → MBDyn proposto

### Matrice richiesta

Usare la base originale `Φ=[φ7 ... φ12]`, con identici segni, unità,
normalizzazione e trasformazioni MPC, e ottenere

\[
G(t)=\Phi^T K_{geo}^{FEM}[\sigma_b(t)]\Phi.
\]

È una matrice 6×6, generalmente non diagonale. Sei frequenze non determinano
i suoi 21 coefficienti indipendenti quando simmetrica. Il confronto di
frequenze fornisce al più una sensibilità diagonale approssimata. Una scala
di frequenza `s` equivale a una scala di rigidezza `s²`, a massa fissa.

Serve esportare la matrice geometrica in un set di DOF noto, non la rigidezza
già alterata dagli SPC del problema modale, e proiettarla nella base compatibile.
Si può equivalentemente sottrarre K precaricata e K scarica solo se set,
vincoli, riduzioni e altri termini coincidono e si identifica cosa contiene
la differenza. Gli attuali output di autovalori e forme su 45 GRID non bastano.

La procedura DMAP/OUTPUT4 va verificata per la versione MSC effettiva su un
caso piccolo e sullo stato già disponibile. Non c'è qui un ALTER universale
già validato. In particolare non si deve chiamare una matrice `KDD` “geometrica”
soltanto perché il nome contiene D: la notazione dei set va verificata.

### Ridurre i calcoli offline usando la baseline

1. Correggere il recupero DLM leggendo i coefficienti dal deck di ogni run;
   controllare risultanti, momenti, reazioni SPC, unità e lavoro virtuale.
2. Recuperare storie di carico, non soltanto le due medie. I NetCDF recovery
   esistenti coprono già la zona di osservazione; non coprono necessariamente
   l'intera preistoria della manovra. Per quest'ultima servirebbe un recupero
   aggiuntivo con output opportuno, non l'intero sweep.
3. Partire da uno stato di riferimento 1 g e tre stati a inizio/centro/fine
   finestra 1.6 g. Sono campioni pilota, non una garanzia di sufficienza.
   Oppure decomporre le storie in pochi pattern di carico indipendenti
   mediante SVD, scalando forze e momenti con una lunghezza dichiarata.
4. Nel regime di statica lineare/stress correction, la matrice geometrica è
   lineare nel campo di sforzo: calcolare le risposte dei pattern indipendenti
   consente di ricombinarle. Riutilizzare la fattorizzazione statica quando
   la sequenza Nastran lo consente, anziché risolvere una modale a ogni tempo.
5. Validare su almeno un carico non usato nella ricostruzione: piccolo errore
   SVD dei carichi non garantisce piccolo errore di G o del polo BFF.
   Raffinare solo se la variazione entro finestra o il residuo lo richiedono.

Il campo DLM fisico resta un'ipotesi separata finché non si recuperano le
pressioni. La sua incertezza va dichiarata, non assorbita in un unico “KG vero”.

## 5. Confronto efficace durante la stessa manovra

Per isolare l'effetto sul BFF propongo prima un esperimento tangente lungo
la shadow originale, non una nuova manovra ri-trimmata. Siano `q_b(t)` la
coordinata della baseline e `G_b(t)` la correzione geometrica ricostruita:

\[
Q_{geo}(t)=-G_b(t)[q(t)-q_b(t)].
\]

La forza è nulla sulla baseline e aggiunge G alla rigidezza delle piccole
perturbazioni. Così l'effetto sulla crescita oscillatoria non si confonde
subito con un cambiamento di trim/comando. Una variazione di G può modificare
sia frequenza sia smorzamento aeroelastico per accoppiamento con il flusso.

**Portata precisa:** è un confronto controllato della dinamica incrementale
condizionato alla traiettoria scelta. Non è il modello nonlineare completo
della struttura precaricata: non ricostruisce il nuovo equilibrio, il feedback
degli sforzi dovuto alla perturbazione o tutti gli accoppiamenti mancanti.
Applicare invece `-G_b(t)q(t)` cambia anche il carico medio; senza un termine
affine coerente non preserva la baseline. Va trattato come esperimento diverso.

L'implementazione può sfruttare forze modali accoppiate, ma deve fornire o
verificare lo Jacobiano della dipendenza da q. Per un uso efficiente e robusto,
un piccolo elemento con residuo `-G(q-q_b)` e derivata `-G` è preferibile a
presumere che generiche drive di forza aggiornino automaticamente lo Jacobiano.
L'interpolazione della baseline deve essere coerente col passo e sufficientemente
accurata: una shadow corretta deve riprodurre la baseline entro l'errore numerico.
Non introdurre rampe/gradini di rigidezza dentro la finestra di identificazione.

### Campagna minima progressiva

| Fase | Nuove run MBDyn | Scopo |
|---|---:|---|
| V=66.75, 1.6 g, dt=0.01, shadow + excited corretti | 2 | Effetto tangente e verifica della baseline |
| V=66.75, 1 g, dt=0.01, stessa coppia | 2 | Separare precarico 1 g ed incremento di manovra |
| Solo la coppia decisiva, dt=0.005 | 2 | Convergenza temporale del confronto |

Totale iniziale **4 run**, poi **2 di verifica** se il segnale è leggibile.
Le coppie originali vengono riutilizzate; non si ripetono le 28 integrazioni
dello sweep. Se necessario, si aggiunge soltanto una verifica di ampiezza
(rap 0.10° invece di 0.20°) o dell'ipotesi di carico DLM sul caso più sensibile.
Il costo in ore non è ancora misurato per il nuovo elemento/AD: promettere lo
stesso runtime delle run correnti sarebbe prematuro. Il risparmio certo è nel
numero di traiettorie e di fattorizzazioni richieste.

Per ciascuno stato confrontare `excited-shadow`, non due segnali assoluti
con trim diverso; misurare Δσ, Δf e la coerenza fra q7 e tip simmetriche.
L'effetto incrementale della manovra è una differenza di differenze:

\[
\Delta\Delta\sigma=
(\sigma_{geo,1.6}-\sigma_{base,1.6})-
(\sigma_{geo,1}-\sigma_{base,1}).
\]

Tenere fissi polari, DLM, M, smorzamento strutturale, SAS, rap e protocollo di
rilascio. Verificare TAS, n, pitch-rate, assetto e superfici; non aumentare
silenziosamente lo smorzamento quando cambia la rigidezza.
La precedente differenza temporale di σ=0.0161 s⁻¹ è un riferimento, non
l'intero intervallo di incertezza del nuovo esperimento. Una finestra di
2.05 s contiene circa quattro cicli: il risultato è un tasso locale medio,
non automaticamente una frontiera di flutter stazionaria.

Solo dopo un effetto superiore a errori di passo, identificazione e mapping,
raffinare la velocità intorno al cambio di segno, senza griglia rettangolare
estesa. Se invece serve la nuova traiettoria fisica, procedere al gruppo 19
con feedback dei carichi e validazione separata del nuovo equilibrio.

## 6. Cosa è stato prodotto in questa revisione

- Questo documento e l'avviso di rettifica nel README.
- `review_existing_results.py`: diagnostica separata, senza cancellazioni o
  modifiche ai deck/NetCDF/F06. Riporta tutte le radici negative senza soglia
  e ricalcola Q7 dai coefficienti dei deck realmente eseguiti.
- `results/model_review/existing_results_audit.json`: risultati numerici con
  percorsi e SHA-256 degli input verificati.

Riproduzione:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 review_existing_results.py
```

Restano da implementare/validare l'esportazione KG, il recupero corretto dei
carichi nei deck futuri e l'elemento/collegamento dinamico. Non sono presentati
come già pronti per simulazioni fisiche. Le eventuali run Nastran rimangono
a carico dell'utente; per nuovi deck/run il percorso concordato è ZENO.
