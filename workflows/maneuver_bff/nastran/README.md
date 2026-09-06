# Gate Nastran per lo stress-stiffening

## Avviso di revisione — 5 settembre 2026

Le conclusioni di stabilità e robustezza del gate riportate di seguito sono
superate dalla [revisione completa](../REVIEW_MODAL_NASTRAN_DLM.md): entrambe
le run a 1.6 g hanno una radice negativa precedentemente filtrata sotto 0.1 Hz,
e il precarico `with_dlm` usa coefficienti diversi dai deck MBDyn eseguiti.
Non usare il precedente gate come validazione fisica o autorizzazione a
trasferire i modi supportati nel FEM free-free. Deck e risultati sono stati
conservati; l'audit separato è `../review_existing_results.py`.

## Decisione attuale dopo lo screening

Le tre campagne MBDyn sono complete. Il gate ha restituito
`physical_stress_stiffening_rom_recommended`: lo span di `sigma` previsto per
una variazione modale da -1% a +1% è 0.22384 1/s a 1 g e 0.24856 1/s a
`n=1.6`, contro una risoluzione di 0.05 1/s. Una nuova elaborazione strutturale
è quindi necessaria se si vuole concludere sull'effetto fisico del prestress,
non per ripetere la campagna lineare già completata.

Il 4 settembre 2026 sono state completate due run MBDyn dedicate di recupero
carichi e sono stati creati quattro deck prestressed fisici: 1 g/1.6 g,
ciascuno con e senza la distribuzione canonica del contributo DLM. Nastran non
è presente nell'ambiente locale, quindi solo queste quattro SOL 103 restano da
eseguire manualmente.

Il caso esistente
`NASTRAN_SIMULATIONS/03_GRAVITY_5G/nastran/MAIN/sol101_gravity_5g.bdf` non va
usato a questo scopo: applica 5 g in -Z e vincola tutti i sei gradi di libertà
del nodo 990001 tramite un RBE2. È un ottimo benchmark statico della base
modale, ma è una semiala/velivolo root-clamped, non un velivolo libero in
equilibrio aero-inerziale.

## Evidenza riproducibile

Dopo le run di `stiffness_screen`, eseguire:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 analyse_sweep.py /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/stiffness_screen \
  --reference-directory /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/primary
```

Il campo `prestress_decision_gate.decision` in
`C:\Users\Utente\Desktop\BFF_PULLUP_V2\stiffness_screen\analysis\summary.json`
vale:

- `prestress_not_resolved_skip_for_primary_conclusion`: l'effetto previsto da
  uno spostamento modale ±1% è inferiore alla risoluzione di 0.05 1/s;
- `physical_stress_stiffening_rom_recommended`: l'effetto è risolvibile e va
  costruito il ROM fisico;
- `insufficient_results`: mancano almeno tre scale valide per una regressione.

Lo screening non stima quale sia lo spostamento reale di frequenza dovuto ai
carichi: misura solo la derivata locale del flutter rispetto alla frequenza del
modo 7.

## Percorso fisico preferito: RECORD GROUP 19

MBDyn descrive la correzione come

\[
K_{qq,geo}=\sum_i(a_{0i}-g_i)K_{0t_i}
 +\sum_i\dot\omega_{0i}K_{0r_i}
 +\sum_m\omega_{q_m}K_{0\omega_m}
 +\sum_{i,j}F_{Pij}K_{0F_{ij}}
 +\sum_{i,j}M_{Pij}K_{0M_{ij}}.
\]

Le matrici ridotte di carico unitario devono essere inserite nel FEM MBDyn nel
`RECORD GROUP 19`. Per la richiamata simmetrica sono prioritarie:

1. le tre `K0t` (tag 10–12), perché la matrice usa l'accelerazione del nodo
   modale meno la gravità;
2. le `K0F` e `K0M` dei 45 nodi di interfaccia realmente caricati (tag da 13
   in poi, sei tag per nodo);
3. `K0r` e `K0omega` solo se le accelerazioni e velocità angolari misurate ne
   rendono significativo il contributo.

Il modello MBDyn deve inoltre contenere `use automatic differentiation;` nel
`control data`: nel sorgente installato solo il ramo `ModalAd` riceve e usa le
matrici del gruppo 19.

La utility `femgen.f90` installata scrive i gruppi 1–11 e non il 19. Occorre
quindi un generatore dedicato, per esempio `mboct-fem-pkg`, oppure una catena
verificata che proietti le matrici geometriche Nastran sulla stessa base modale
e scriva il formato MBDyn. Non basta rieseguire l'attuale SOL 103 e passarlo a
`femgen`.

## Verifica disponibile: snapshot prestressed SOL 103

Le shadow di recupero sono a 66.75 m/s, passo 0.01 s, senza rap. Sono mediate
nelle finestre interne SAS-off 13.25–14.79 s a 1 g e 16.38–17.90 s a 1.6 g.
Il recupero integra i tre punti di Gauss dei 58 corpi aerodinamici e trasporta
forze e momenti sui 45 reference grid RBE3. Il fattore di carico medio è 0.9940
e 1.6329; la portanza recuperata è 417.278 e 694.622 lbf.

L'inerzia traslazionale è applicata con `GRAV=M(g-a_CG)` e quella rotazionale
con `RFORCE`, includendo sia velocità sia accelerazione angolare. Il vincolo
`SPC1,620,123456,10062` è selezionato sia nel subcase statico sia nel subcase
modale con `STATSUB=1`. Si estraggono così i modi elastici attorno allo stato
caricato con un supporto comune al CG. Questa scelta separa lo stress stiffening
strutturale dalle radici rigide spurie generate dai carichi esterni dead in una
linearizzazione free-free; segue inoltre l'esempio SOL 103 prestiffened della
guida MSC, che mantiene lo stesso SPC nei due subcase. Una forza/coppia sullo
stesso grid vincolato chiude le risultanti e
rappresenta la reazione X/Y del pin MBDyn. Il residuo è quasi tutto assiale
(33.22 lbf a 1 g, 90.31 lbf a 1.6 g); quello verticale è 0.364/9.718 lbf e il
residuo in momento ha norma 0.285/1.003 lbf-in.

Ogni cartella contiene:

- `preload_loads.bdf`: carichi aero, `GRAV`, `RFORCE`, bilanciamento e DLM;
- `prestressed_modes_nodlm.bdf`: verifica senza distribuzione DLM;
- `prestressed_modes_with_dlm.bdf`: verifica con
  `M_lumped*phi_7*Q_7`, normalizzata a proiezione modale unitaria;
- `load_summary.json` e `nodal_loads.csv`: tracciabilità numerica;
- `run_nastran.sh`: lancio sequenziale dei due deck.

La cartella `n1p0` contiene anche `supported_unloaded.bdf`: è l'unica run di
riferimento aggiuntiva, senza carichi e con lo stesso `SPC=620`. Serve a
separare la variazione dovuta al precarico dalla variazione apparente introdotta
dal passaggio tra modi baseline free-free e modi supportati al CG.

La variante DLM è un bracket, non un recupero di carichi di pannello: la forza
generalizzata media è +0.7951 a 1 g e -60.1338 a 1.6 g e non determina una
distribuzione fisica unica.

Esecuzione:

```bash
cd /mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2/load_recovery/prestress_loads/n1p0
./run_nastran.sh
cd ../n1p6
./run_nastran.sh
```

Il wrapper usa `${NASTRAN_CMD:-nastran}` e presuppone l'interfaccia standard
`nastran deck.bdf scr=yes old=no`. Sul server ZENO il wrapper `nast` usa invece
l'ordine `nast old=no deck.bdf`: in quel caso lanciare manualmente, aspettando
`END OF JOB` prima del comando seguente:

```bash
cd ~/ZENO/prestress_loads
./run_all_zeno.sh
```

Il generatore omette separatamente le card `FORCE` e `MOMENT` con direzione
nulla. Questo controllo evita `USER FATAL 9994 (BULFUN)` nei nodi in cui la
realizzazione modale DLM possiede soltanto componenti traslazionali oppure
soltanto componenti rotazionali.

I deck specificano inoltre `PARAM,TESTNEG,3`. Per la Quick Reference MSC questo
valore mantiene la differential stiffness anche se la decomposizione incontra
termini negativi, consentendo l'estrazione degli autovalori anziché la chiusura
con `UFM 4413`. Non è una soppressione del controllo fisico: il postprocess
elenca tutte le radici negative e, se ne trova anche una nel problema vincolato
al CG, classifica lo stato come instabile/non valido e impedisce la costruzione
del ROM.

Controllare in ciascun `.f06`:

- assenza di `USER FATAL` e presenza di `END OF JOB`;
- completamento sia del subcase statico sia di quello modale;
- assenza di modi rigidi nel subcase 2 supportato al CG;
- autovalori e autovettori stampati sui 45 grid del `SET 901`.

Poi eseguire:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 analyse_prestressed_modes.py
```

Il codice abbina i modi 7–12 tramite MAC. Il ROM viene raccomandato soltanto se
le due varianti hanno `MAC>=0.90`, uno shift 1 g→1.6 g del modo 7 almeno pari a
0.4% e lo stesso segno. Un disaccordo tra i bracket impone di recuperare i
carichi DLM di pannello prima di proseguire.

Questa verifica congela il precarico a un solo stato e non aggiorna `Kgeo(t)`
durante la richiamata. È una validazione frozen-time, non un sostituto completo
del gruppo 19. Il thrust assiale resta una reazione al CG perché questa è la
fisica del modello MBDyn usato nella campagna.

## Risultati finali

Tutte le cinque run terminano con `END OF JOB` e non presentano radici elastiche
negative. Il modo FEM 7 è tracciato sul modo supportato 4 con MAC circa 0.920.
Con riferimento scarico supportato di 2.768387 Hz, a 1.6 g si ottengono
2.752805 Hz senza DLM (-0.5629%) e 2.751967 Hz con DLM (-0.5931%). A 1 g gli
shift sono rispettivamente -0.1850% e -0.1847%. Il bracket concorda quindi su
un softening superiore alla soglia 0.4% nello stato di richiamata e il gate
richiede un ROM fisicamente stress-stiffened. I modi supportati di questo
screening non sostituiscono direttamente la base free-free MBDyn: il passo
successivo corretto è la costruzione delle matrici ridotte del gruppo 19.

Riferimenti ufficiali usati per il percorso Nastran:

- [MSC Nastran Dynamic Analysis User's Guide 2022.2](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.2_Dynamic_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.2_Dynamic_Analysis_User_Guide.pdf), sezione sui modi preloaded con `STATSUB(PRELOAD)`;
- [MSC Nastran Aeroelastic Analysis User's Guide 2022.3](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.3_Aeroelastic_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.3_Aeroelastic_Analysis_User_Guide.pdf?save_local=true), per trim elastico e recupero dei carichi.
