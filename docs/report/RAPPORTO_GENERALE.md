# Rapporto generale della cartella TESI

**Snapshot:** 5 settembre 2026, Europe/Rome  
**Radice:** `/home/nicomonzi/TESI`  
**Oggetto:** modello aeroelastico NASA X-56/MUTT, validazioni Nastran–MBDyn,
trim e controllo, identificazione del body-freedom flutter (BFF), manovra e
accoppiamento MBDyn–DUST.

## Sintesi esecutiva

Il progetto ha superato la fase di puro prototipo: esistono baseline
strutturali e aerodinamiche, workflow riproducibili, controlli di qualità,
campagne concluse e confronti numerici. Non è però ancora corretto definirlo
“concluso” o “validato in produzione”. Le tre linee più mature sono:

1. `BFF_open_loop`, baseline MBDyn open-loop con SAS-off, che ha mostrato un
   cambio stabile/instabile tra 65 e 70 m/s e una stima preliminare di frontiera
   pari a circa 66.34 m/s;
2. `NASTRAN_SIMULATIONS`, che ha costruito la base a 60 modi, misurato la
   sensibilità COUPMASS e completato lo sweep statico 5 g da 1 a 54 modi;
3. `MANOUVER_STIFNESS`, il ramo più recente, che ha eseguito 28 traiettorie più
   due recuperi di carico e cinque SOL 103 lato Zeno, ma ha poi rettificato le
   proprie conclusioni: radici negative e incoerenza dei coefficienti DLM
   rendono non validato il gate fisico sul prestress; inoltre il nuovo stimatore
   non conferma ancora i vecchi onset di manovra senza i controlli causali.

L'accoppiamento `BFF_DUST_55` è tecnicamente ben organizzato in CHECK/SMOKE/
PRODUCTION. Lo smoke da 50 passi è completato; la produzione resta bloccata da
mesh senza winglet, mancata convergenza formale, hinge non allineate, conflitto
sulla base modale e verifica DUST dei controlli non eseguita.

Il rischio operativo principale non è il codice: è la gestione degli asset.
La radice contiene circa 5.15 GB di file di lavoro, mentre `.git` occupa circa
5.82 GiB in oggetti. Sono presenti ambienti virtuali, build, cache, output grezzi
e grandi FEM replicati. Il ramo scientifico più recente,
`MANOUVER_STIFNESS`, è interamente non tracciato da Git nello snapshot.

## Numeri dello snapshot

| Indicatore | Valore |
|---|---:|
| File regolari analizzati | 13.256 |
| Dimensione dei file analizzati | 5.154.403.381 byte, circa 4,80 GiB |
| Link simbolici | 4 |
| File vuoti | 252 |
| Commit Git | 68 |
| Periodo della cronologia Git | 7 maggio–28 agosto 2026 |
| Oggetti Git sciolti | 5.415, 1,90 GiB |
| Oggetti Git in pack | 16.910, 3,92 GiB |
| Garbage Git | 6 file, 35,71 MiB |
| Gruppi di file duplicati non vuoti | 383 |
| File coinvolti nei gruppi | 2.427 |
| Copie oltre la prima | **2.044** |
| Recupero massimo teorico | 1.132.921.243 byte, circa 1,06 GiB |

I dati completi sono in `dati/inventario_file.csv`,
`dati/duplicati_sha256.csv`, `dati/riepilogo_cartelle.csv`,
`dati/inventario_sottocartelle.csv` e `dati/catalogo_codice_python.csv`. I 17
rapporti interpretativi corrispondono alle sottocartelle di primo livello;
l'inventario ricorsivo copre anche ogni directory annidata.

## Architettura scientifica ricostruita

```text
Nastran BDF/OP2/F06
    │  SOL 103, SOL 101, SOL 145; femgen
    ▼
FEM modale MBDyn ──► modal joint floating-frame ──► modello strutturale X-56
    │                                                     │
    │                                   C81/Theodorsen + correzione DLM
    │                                                     │
    ├──► trim 2×2 Fz/My ──► SAS/attuatori/filtri ──► RAP + SAS-off
    │                                                     │
    │                                  Matrix Pencil / fit oscillatorio
    │                                                     ▼
    ├──► sweep BFF steady                         onset e tasso di crescita
    └──► coppie shadow/excited durante dive–pull-up ──► effetto della manovra

MBDyn nodal socket ◄──── preCICE ────► DUST 3-D panels/wake
                         │
                         └──► CHECK, smoke, diagnostica e futura produzione
```

La scelta metodologica ricorrente è buona: separare generazione input,
esecuzione e post-processing; conservare l'input renderizzato; rifiutare run
incomplete; usare manifest e SHA-256; confrontare segnali indipendenti; non
trasformare uno smoke o uno screening parametrico in una conclusione fisica.

## Stato cartella per cartella

| Cartella | File / dimensione apparente | Ruolo | Stato | Obsoleta? |
|---|---:|---|---|---|
| `BFF_DUST_55` | 1.750 / 229 MB | Accoppiamento MBDyn–DUST finale | Attiva, produzione bloccata | No; output visuali in parte rigenerabili |
| `BFF_maneuver_envelope` | 21 / 148 MB | Generatore dive–pull-up e analisi paired | Attiva, modifiche non consolidate | No; è dipendenza del ramo manovra |
| `BFF_open_loop` | 35 / 56 MB | Baseline BFF steady SAS-off | Canonica e più matura | No |
| `DUST` | 15 / 21 kB | Geometria/airfoil e input DUST originari | Supporto corrente | No, ma va centralizzata |
| `GRAFICI` | 349 / 150 MB | Figure di validazione e materiale tesi | Deliverable/supporto | No; separare sorgenti e artefatti |
| `MANOUVER_STIFNESS` | 195 / 5,5 MB | Studio più recente di manovra/prestress | Attiva con conclusioni rettificate | No; `runs_superseded_full_matrix` sì |
| `MBDYN` | 217 / 784 MB | Sviluppo storico: massa, SOL101, PID, autopilot, V-g | Archivio misto | In gran parte sì per esecuzione corrente |
| `NASTRAN` | 8.885 / 559 MB | Modelli originali/rilasciati, parser e viewer | Riferimento misto | Non tutta; `.venv_x56` e cache sì |
| `NASTRAN_SIMULATIONS` | 643 / 455 MB | Validazioni modali, COUPMASS e 5 g | Supporto attivo, validazione esterna aperta | No |
| `TEST` | 846 / 1,8 GB | Prototipi trim/BFF/DUST e relativi output | Archivio sperimentale | Sì come linea attiva; preservare evidenze |
| `TRIM` | 43 / 189 MB | Trim longitudinale MIMO 2×2 e Jacobiano | Sottosistema validato | No; NetCDF di calibrazione archiviabili |
| `X56_AERO_POLAR` | 131 / 126 MB | Confronto polar Nastran–MBDyn e C81 spanwise | Studio storico/supporto | Parzialmente superata dalla C81 effettiva |
| `_trim_coupled_stage` | 31 / 55 MB | Staging trim MBDyn–DUST a 55 m/s | Preparata, senza risultato finale locale | Non ancora; incubatore da integrare |
| `bbf_manouver` | 44 / 245 MB | Primo workflow BFF durante manovra | Precursore con casi validati | Sì, sostituito dai rami dive–pull-up |
| `bff_eigenalaisy` | 28 / 64 MB | Studio 4-DOF closed-loop | Un caso inconclusivo | Sì, sostituito |
| `bff_longitudinal` | 18 / 56 MB | SAS trasparente con notch multipli | Preparato ma non simulato | Sì, sostituito da `BFF_open_loop` |
| `results` | 8 / 1,0 MB | Grafici/CSV globali di confronto | Output fuori sede | Non scientificamente; sì come collocazione |

Le cartelle nascoste `.agents` e `.codex` sono vuote. I rapporti dettagliati
sono in `cartelle/`.

## Cosa è stato fatto e come

### 1. Fondazione strutturale Nastran

Tra maggio e luglio sono stati importati e corretti i deck X-56, prodotti
risultati SOL 101/103/145 e generati file FEM per MBDyn. La catena usa BDF/BULK
come sorgente, OP2/MAT/F06 come output e `femgen` per esportare nodi, forme,
masse e rigidezze modali. Sono nati anche parser e viewer per verificare forme,
frequenze, massa e deformate.

La validazione più ordinata è ora in `NASTRAN_SIMULATIONS`: 60 modi totali,
sei rigidi e 54 elastici; 45 coppie COUPMASS affidabili su 54; convergenza
interna MBDyn raggiunta con 30 modi elastici come scelta robusta. Rimane uno
scarto esterno alle tip di circa 1,17 in con 54 modi. La causa più plausibile
documentata è l'incoerenza RBE3/RBE2 tra FEM modale e SOL 101, aggravata dalla
grande deformazione lineare a 5 g. Quindi la convergenza interna è dimostrata,
la validazione Nastran–MBDyn all'1% no.

### 2. Controllo, trim e attuatori MBDyn

A giugno sono stati sviluppati PID di rollio e beccheggio, autopilota, filtri
digitali, notch e attuatore del primo ordine. Il ramo `TRIM` ha poi formalizzato
un trim longitudinale MIMO con residui `[Fz, My]` e incognite `[pitch,
body-flap]`. Il Jacobiano 2×2 è stato identificato con differenze finite
centrate; il passo Newton di validazione ha ridotto i residui da
`[-71,10 lbf, 631,29 lbf in]` a `[1,52 lbf, 1,10 lbf in]`.

Il modello finale separa anelli lenti e dinamica aeroelastica: quota/Vz sulle
superfici alari, pitch/q sui body flap, roll/p in differenziale e filtri
discretizzati in modo coerente col passo. Questa separazione ha corretto il
doppio percorso quota–pitch che aveva prodotto un ciclo limite e divergenza.

### 3. Identificazione BFF steady

I prototipi `TEST/bbf`, `TEST/BBF_AUTO_TRIM`, `bff_eigenalaisy` e
`bff_longitudinal` hanno progressivamente introdotto volo libero 4-DOF,
controllo, burst WF4, analisi spettrale/Hilbert/Matrix Pencil e trasparenza del
SAS. La sintesi corrente è `BFF_open_loop`:

- ingresso controllato al punto;
- controllo automatico dello stato prima del rilascio;
- hold di tutte le superfici e azzeramento reale del feedback per 2,05 s;
- piccolo RAP simmetrico WF4;
- identificazione congiunta di tip, modo FEM 7 e velocità modale;
- riattivazione del SAS e recupero.

I casi documentati a 57,5 e 65 m/s sono stabili, quello a 70 m/s instabile. La
frontiera lineare preliminare è circa 66,34 m/s, vicina al riferimento SOL 145
di 65,764 m/s TAS. La frequenza MBDyn a 70 m/s è però sottostimata, segnale che
la correzione DLM monomodale non è quantitativamente completa.

### 4. Aerodinamica 3-D DUST e accoppiamento

Da luglio la geometria parametrica DUST è stata dotata di superfici mobili e
collegata a MBDyn tramite socket nodale e preCICE. Gli adapter scambiano
posizione/rotazione/velocità e carichi; i runner generano geometria, validano
gli eseguibili, avviano i due solver, controllano la terminazione e preparano
VTU/PVD per ParaView.

Il ramo finale `BFF_DUST_55` conserva configurazioni preCICE 2.x e 3.x, tre
livelli di mesh e una macchina a stati CAPTURE_TRIM → READY → RAP → OPEN_LOOP →
RECOVERY. Lo smoke COARSE da 900 pannelli, 59 nodi, `dt=0,002 s` e 50 finestre
è terminato in circa 127 s. La FINE da 4.284 pannelli è solo una mesh
diagnostica: mancano winglet, studio originale e allineamento hinge; le
differenze MEDIUM→FINE restano 3,16% in Fz, 14,59% in My e 26,33% nella
distribuzione spanwise.

### 5. BFF durante manovra e prestress

`bbf_manouver` ha verificato casi level, pull-up e roll con SAS sempre attivo.
`BFF_maneuver_envelope` ha poi introdotto coppie `shadow/excited`, una vera
finestra SAS-off e una dive–pull-up parametrica. `MANOUVER_STIFNESS` ha ridotto
una matrice iniziale di 78 traiettorie a 28 casi mirati, con riduzione del 64%
dei casi e del 79% dei passi:

- 18 traiettorie principali su 66,25/66,75/67,25 m/s e n=1,0/1,3/1,6;
- 2 per convergenza temporale 0,01/0,005 s;
- 8 per sensibilità della frequenza del modo 7 di ±1%;
- 2 shadow aggiuntive per recupero carichi 1 g e 1,6 g.

Il confronto paired sottrae `excited-shadow`, così la manovra comune viene
rimossa al primo ordine. Il vecchio stimatore Hilbert indicava onset
66,902/66,834/66,810 m/s e convergenza temporale; lo screening mostrava forte
sensibilità a K77. La revisione successiva ha però mostrato che il fit Hilbert
può invertire o distorcere il segno su finestre corte. Il fit diretto
polinomio+sinusoide smorzata trova risultati discordanti e richiede tre casi
causali aggiuntivi. Pertanto quei tre onset sono risultati storici, non la
conclusione corrente.

Per il prestress, le SOL 103 frozen-time mostrano un softening incrementale del
modo 7 circa -0,379% senza DLM e -0,409% con DLM, ma compaiono radici negative a
1,6 g e il recupero DLM usa coefficienti incoerenti. Il modal joint non contiene
il RECORD GROUP 19, quindi non implementa stress-stiffening dinamico. La strada
corretta è diagnosticare le radici, correggere il recupero e costruire una
matrice geometrica ridotta nella stessa base free-free, non inserire modi
supportati direttamente nel modal joint.

## Duplicati: interpretazione e priorità

I due gruppi dominanti sono versioni diverse della base FEM modale:

- `mbdyn_modal.fem`, 38.715.261 byte, 15 copie: massimo teorico circa 516,9 MiB;
- `mbdyn_modal_60.fem`, 57.377.461 byte, 8 copie: massimo teorico circa
  383,0 MiB;
- `nsvibe_test.fem`, 10.775.961 byte, 14 copie: massimo teorico circa 133,6 MiB.

Questi tre gruppi spiegano quasi tutto il recupero teorico. Sono copie
intenzionali usate per rendere i casi autonomi; eliminarle senza cambiare gli
include romperebbe le simulazioni. La soluzione non è cancellare alla cieca,
ma creare un archivio di asset versionati per hash e far puntare i workflow a
una versione esplicita.

Duplicati più sicuri da eliminare o rigenerare dopo backup:

- `__pycache__` e 3.605 file `.pyc`, circa 82 MB;
- `.venv_x56` dentro `NASTRAN`, circa 348 MB;
- build locale `TEST/BBF_DUST_COUPLED/build_dust`, circa 33 MB;
- sei `incomplete_attempt_*` del coupling test;
- frame neutri e file airfoil copiati in decine di directory di animazione;
- script chiamati `analyse_bff copy.py`, identici in tre rami;
- output `.nc`, `.aer`, `.log`, `.out`, `.vtu` quando esiste input, versione
  solver e manifest sufficiente a rigenerarli.

Non vanno eliminati automaticamente: F06/OP2/MAT di riferimento, NetCDF che
costituiscono l'unica evidenza di una run costosa, input renderizzati usati per
una conclusione, manifest e CSV/JSON finali. Il catalogo CSV permette una
revisione gruppo per gruppo.

## Cartelle obsolete o superseded

### Obsolescenza esplicita

- `MANOUVER_STIFNESS/runs_superseded_full_matrix`: il README dice di non
  eseguirla; è la matrice estesa sostituita da quella a 28 traiettorie.
- ogni `__pycache__`, `.pyc`, `.venv_x56` e `build_dust`: artefatti di ambiente,
  non sorgenti scientifiche.
- `TEST/BBF_DUST_COUPLED/output/.../incomplete_attempt_*`: tentativi falliti
  conservati automaticamente.

### Rami scientifici superati come esecuzione corrente

- `bff_eigenalaisy` → `bff_longitudinal` → `BFF_open_loop`;
- `bbf_manouver` → `BFF_maneuver_envelope` → `MANOUVER_STIFNESS`;
- `TEST/bbf` e `TEST/BBF_AUTO_TRIM` → workflow BFF successivi;
- `TEST/BBF_DUST_COUPLED` → `_trim_coupled_stage`/`BFF_DUST_55`;
- `MBDYN/PID*`, `AUTOPILOT`, `BBF` → controlli integrati nelle baseline
  successive;
- parte di `X56_AERO_POLAR/test` → `BFF_open_loop/build_x56_c81.py` e relative
  validazioni statiche.

Questi rami vanno archiviati, non cancellati: contengono la genealogia delle
scelte e talvolta l'unica run disponibile. Prima dell'archiviazione occorre
aggiungere un file `STATUS.md` e una mappa “sostituito da”.

## Riorganizzazione proposta

```text
TESI/
├── docs/
│   ├── report/                 # questo dossier
│   ├── metodologia/
│   └── decision_log/
├── models/
│   ├── nastran/x56_r11/
│   ├── mbdyn/x56_modal/
│   ├── aerodynamics/c81/
│   └── dust/x56/
├── workflows/
│   ├── trim/
│   ├── bff_open_loop/
│   ├── maneuver_bff/
│   └── coupled_dust/
├── validation/
│   ├── modal/
│   ├── static_5g/
│   ├── aero_polar/
│   └── coupling_mesh/
├── results/
│   ├── manifests/
│   ├── processed/
│   └── figures/
├── assets/                     # FEM/BULK grandi, versionati per SHA-256
├── tools/
└── archive/
    ├── legacy_mbdyn/
    ├── prototypes_test/
    └── superseded_campaigns/
```

La migrazione dovrebbe avvenire in quattro fasi, sempre con commit separati:

1. congelare lo stato corrente, aggiungere `MANOUVER_STIFNESS` a Git e salvare
   un manifest degli output esterni ZENO/Desktop;
2. introdurre un registro canonico degli asset grandi con nome, SHA-256,
   origine, unità e workflow consumatori, senza spostare ancora nulla;
3. spostare prima cache/build/output rigenerabili, poi codice, aggiornando e
   verificando include e import dopo ogni gruppo;
4. archiviare i rami superati con README immutabile e mantenere nei workflow
   attivi solo link/configurazioni verso asset canonici.

## Priorità operative

1. **Mettere in sicurezza il lavoro corrente:** versionare o almeno fare backup
   di `MANOUVER_STIFNESS` e delle modifiche in `BFF_maneuver_envelope`.
2. **Chiudere il test causale:** completare i tre controlli descritti in
   `SWEEP_PULLUP_BFF.md` e rigenerare l'analisi col fit diretto.
3. **Correggere il prestress:** diagnosticare le radici negative, allineare i
   coefficienti DLM e decidere se costruire RECORD GROUP 19/Kgeo ridotta.
4. **Chiudere la validazione strutturale:** rigenerare il FEM con condizioni
   RBE2 coerenti col SOL 101 e ripetere il confronto a 5 g o a carico lineare
   più basso.
5. **Sbloccare DUST solo dopo i gate:** winglet/topologia, hinge, convergenza,
   base modale ed efficacia controlli; poi breve run server per stimare memoria
   e particelle prima delle 4.750 finestre.
6. **Ridurre il debito di struttura:** eliminare solo cache/venv/build dopo
   backup, quindi centralizzare FEM e BULK senza perdere riproducibilità.

## Aggiornamento continuo

`report/aggiorna_report.py` rigenera inventario, duplicati, riepilogo e
cronologia Git. Confronta gli SHA-256 con lo snapshot precedente e aggiunge a
`PROGRESSI.md` file aggiunti, rimossi e modificati. I rapporti interpretativi
non vengono riscritti automaticamente: quando cambia una conclusione
scientifica, occorre aggiornare questo file e il rapporto della cartella
interessata, citando il nuovo commit o manifest.
