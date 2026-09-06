# Rapporto: `TRIM`

## Ruolo e stato

È il sottosistema dedicato al trim longitudinale MBDyn. Contiene 43 file e
circa 189 MB, quasi tutti dovuti ai sei NetCDF di calibrazione del Jacobiano.
È **attivo come riferimento validato**, anche se i workflow BFF incorporano poi
il controllo nel proprio modello.

## Codice e modello

- `main_trim.mbd` e `INCLUDE/`: modello, nodi, joint e aerodinamica;
- `run_trim_sweep.py`: griglia esatta 30–70 m/s, rilevamento dell'eseguibile,
  render di V_INF/RHO, riuso dei casi completi e manifest incrementale;
- `analyze_trim_sweep.py`: lettura NetCDF, conversioni SI, convergenza, pendenze,
  CSV e figure tesi;
- `identify_trim_jacobian.py`: baseline, perturbazioni centrate pitch ±0,25° e
  body-flap ±0,50°, media finale e validazione Newton;
- `jacobian_calibration_bfl_bfr/`: sei casi, CSV e JSON riproducibile.

Il controllore non è un autopilota classico: risolve dinamicamente
`[Fz, My]=0` usando `[pitch, BFL/BFR]`. L'inversa del Jacobiano disaccoppia i
due canali; i PID hanno solo azione integrale con risposta 0,12 1/s. Un
Butterworth del secondo ordine a 0,25 Hz filtra le reazioni flessibili prima
dell'integratore, e una rampa evita l'inserimento impulsivo.

## Risultati

Il Jacobiano identificato a 63 m/s ha condition number 51,62. Il passo Newton
ha portato il candidato a pitch 0,2488° e body flap 3,1539°. I residui sono
passati da `[-71,10 lbf, 631,29 lbf in]` a `[1,52 lbf, 1,10 lbf in]`, inferiori
alla deviazione standard finale riportata. I coefficienti sono stati inseriti
nel modello principale.

Il README prepara una sweep di 17 casi fra 30 e 70 m/s con output sul Desktop;
quei risultati esterni non sono presenti in questa cartella e quindi non sono
certificati dall'inventario locale.

## Obsolescenza e azioni

Non è obsoleta. Le sei storie NetCDF da 32,75 MB ciascuna possono essere
spostate in un archivio risultati dopo aver conservato JSON, input esatti,
versione MBDyn e hash. Il codice dovrebbe diventare una libreria/validazione del
trim usata dai workflow, evitando implementazioni duplicate nei rami TEST.

