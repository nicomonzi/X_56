# Rapporto: `MANOUVER_STIFNESS`

## Ruolo e stato

È il ramo scientifico più recente: studia l'effetto della dive–pull-up
sull'onset BFF e la sensibilità allo stress-stiffening. Ha 195 file e circa 5,5
MB. È **attivo con blocchi critici** ed è interamente non tracciato da Git nello
snapshot; questa è la priorità di conservazione numero uno.

## Struttura e codice

- `campaign_config.json`: griglia, soglie e parametri unici;
- `campaign.py`: rendering dei casi, ridiscretizzazione coerente di filtri e
  attuatori, variazione parametrica K77 e fingerprint delle sorgenti;
- `run_sweep.py`: prepara/esegue `primary`, `timestep` e `stiffness_screen`;
- `analyse_sweep.py`: paired analysis, qualità, onset e sensibilità;
- `audit_modal_joint.py` e `verify_setup.py`: audit FEM/RECORD GROUP e preflight;
- `run_load_recovery.py` e `extract_prestress_loads.py`: due shadow ad output
  esteso, trasformazione body-frame, integrazione Gauss ed emissione deck SOL
  103 con/senza distribuzione DLM;
- `analyse_prestressed_modes.py`: matching MAC sui 45 nodi di interfaccia;
- `review_existing_results.py`: revisione read-only dei risultati Nastran;
- `run_pullup_causal_controls.py` e `analyse_pullup_bff.py`: nuova verifica
  causale e fit diretto polinomio+sinusoide;
- `runs/`: 28 input ottimizzati più due load recovery;
- `runs_superseded_full_matrix/`: matrice da 78 traiettorie, esplicitamente da
  non eseguire;
- `sweep_pullup_v3/`: nuova estensione preparata con controlli causali;
- `results/` e `audit/`: risultati sintetici e gate.

## Lavoro completato

Sono terminate 18 run primary, due timestep e otto stiffness screen, più due
shadow di recupero: nessun fallimento e copertura paired completa. La campagna
ottimizzata ha ridotto i casi del 64% e i passi del 79%. Il vecchio metodo ha
prodotto onset interpolati 66,9022/66,8341/66,8098 m/s per n=1/1,3/1,6 e ha
verificato `dt=0,01 s` contro 0,005 s con differenza sigma 0,0161 1/s.

Lo screening ±1% del modo 7 ha mostrato span di sigma 0,224–0,249 1/s, quindi
una potenziale sensibilità fisicamente importante. Sono stati recuperati
carichi medi a 1 g e 1,6 g, generati deck e svolte cinque SOL 103 su Zeno.

## Rettifica corrente

Le conclusioni precedenti non sono più il verdetto finale. Il fit diretto ha
dimostrato con segnali sintetici che lo stimatore Hilbert su finestra corta può
essere fortemente dipendente dalla fase. I vecchi onset non sono considerati
verificati; tre controlli causali `sham_release`, `pullup_sas_continuous` e
`sham_sas_continuous` sono ancora necessari.

Il postprocess prestress precedente scartava radici equivalenti sotto 0,1 Hz.
Le run a 1,6 g contengono invece autovalori negativi `-0,0393` e `-0,0796
s^-2`. Inoltre la ricostruzione DLM usa coefficienti diversi dai deck MBDyn.
Lo shift del modo 7 (-0,379%/-0,409%) è quindi indicativo, non un gate validato.
Il FEM contiene RECORD GROUP 1–11 ma non 19: il modal joint usa rigidezza
lineare costante e non stress-stiffening dinamico.

## Obsolescenza e azioni

Solo `runs_superseded_full_matrix` è certamente obsoleta per l'esecuzione.
Conservare i vecchi JSON come evidenza della rettifica, ma marcarli
`superseded`: il campo `decision` nel JSON prestress non riflette la revisione
più recente. Versionare l'intera cartella, completare i controlli causali,
correggere il recupero DLM, diagnosticare le forme negative e costruire Kgeo
nella stessa base free-free prima di qualunque conclusione fisica.

