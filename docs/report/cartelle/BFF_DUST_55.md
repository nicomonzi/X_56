# Rapporto: `BFF_DUST_55`

## Ruolo e stato

È il workflow MBDyn–DUST più recente e meglio impacchettato. Organizza un caso
X-56 a 60,8421 m/s in tre livelli: CHECK, SMOKE e PRODUCTION. È **attivo**, ma
la produzione è intenzionalmente bloccata dal runner; lo smoke disponibile non
è fisicamente trimmato e non può sostenere conclusioni sul flutter.

Nello snapshot contiene 1.750 file e circa 229 MB. Solo 48 file sono tracciati
da Git: i restanti 1.702 sono soprattutto output, cache e prodotti visuali.

## Struttura

- `model/`: modello di produzione, configurazioni preCICE 2.x/3.x, input DUST e
  MBDyn e copia del FEM modale;
- `config/`: esempio di configurazione macchina; separa i path degli eseguibili
  dalla fisica del caso;
- `meshes/`: tre discretizzazioni COARSE/MEDIUM/FINE ricostruite;
- `adapters/`: ponte Python fra socket nodale MBDyn e API preCICE;
- `tools/`: audit mesh, configurazione discretizzazioni, estrazione dei testi
  MBDyn, correzione VTU e analisi di efficacia dei controlli;
- `reference/`: base modale, macchina a stati e piano dei casi di efficacia;
- `reports/`: audit numerico e VTU della mesh FINE;
- `output/`: smoke e check già generati;
- `CASE_DUST_SMOKE/`: caso qualitativo storico a 1.296 pannelli e animazioni
  statiche delle superfici;
- `CASE_DUST_FLUTTER_EMULATION/`: singolo caso open-loop a 3 s, non una sweep;
- `precice-run/`: directory runtime vuota nello snapshot.

## Codice e funzionamento

`run_case.py` è l'entry point. Legge `machine.env`, risolve eseguibili, riconosce
la major version di preCICE, imposta thread/NUMA, genera gli input e applica i
gate. CHECK verifica dipendenze, percorsi, COARSE, nodi di coupling e parsing dei
due solver senza avanzare nel tempo. SMOKE lancia entrambi i processi, controlla
la terminazione e post-processa. PRODUCTION usa una macchina a stati che
mantiene DUST e la scia continui durante trim, RAP, open-loop e recupero.

`coupling_adapter.py` scambia cinematica nodale e carichi; supporta checkpoint e
iterazioni implicite. `extract_mbdyn_text.py` usa `.mov/.mod/.jnt` perché il
NetCDF può restare vuoto durante il coupling iterativo. `audit_mesh.py` calcola
area, aspect ratio, normali, simmetria, bordi non-manifold e allineamento hinge,
poi scrive JSON e VTU.

## Risultati e limiti

Lo smoke finale disponibile usa 900 pannelli, 1.001 punti, 59 nodi, `dt=0,002
s`, 50 finestre e otto thread; è terminato in 127,0 s. La mesh FINE ha 4.284
pannelli e 4.445 punti, senza duplicati o bordi non-manifold, ma solo il 98,60%
delle normali supera il test locale.

I blocchi sono sei: assenza winglet; ricostruzione senza input originali dello
studio mesh; hinge non tutte allineate; differenze MEDIUM→FINE ancora elevate;
conflitto fra numero di modi e frequenza massima; diagnostica rigida di
efficacia DUST non eseguita. Anche preCICE 3.x deve essere validato sul server.

## Obsolescenza e azioni

La cartella non è obsoleta. Sono rigenerabili `__pycache__`, gran parte di
`output/check*`, i frame delle animazioni e molte copie di airfoil/VTU. Le copie
locali del FEM e degli include sono invece dipendenze reali finché non esiste un
asset store canonico. Prima di una produzione: risolvere tutti i gate, misurare
su server il rateo HDF5/particelle con una run breve e solo allora autorizzare le
4.750 finestre previste.

