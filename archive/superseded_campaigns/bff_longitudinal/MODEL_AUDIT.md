# Audit del modello

- Nodo modale: `990000`; giunto modale: `5`. Il total-pin `1` ha due componenti attive: traslazioni globali X e Y.
- Gradi rigidi liberi: Z, roll, pitch e yaw. X e Y sono le sole componenti traslazionali `active`; Z e tutte le rotazioni del joint `1` sono `inactive`.
- Modi FEM attivi: 7–31, cioè i primi 25 elastici. Modi 1–6: rigidi esclusi. Campo di frequenza secco: 3.2171–36.9701 Hz. Il FEM locale contiene 60 modi totali.
- Superfici: BFL/BFR `1004/2004`; WF1 `1008/2008`; WF2 `1011/2011`; WF3 `1014/2014`; WF4 `1017/2017`.
- Combinazione simmetrica: stesso segno numerico L/R. Combinazione differenziale: segno opposto L/R.
- Rappresentazione del moto: velocità strutturale iniziale nulla e vento uniforme `+VINF`; `Vrel` viene sempre ricalcolata come vento meno velocità inerziale.
- Controllo: quota/Vz → riferimento pitch; pitch/q → body flap simmetrici; roll/p → WF1/WF2 differenziali; yaw/r/Vy → body flap differenziali.
- Non esistono propulsione, throttle, airspeed hold o feedback delle reazioni X/Y.
- Il trim e' un termine separato; le correzioni degli attuatori partono da zero e non esiste alcun rilascio programmato del SAS.
- Burst: WF4 simmetrico, finestra Hann, 4 cicli a 2.0597 Hz. Il SAS resta attivo; pitch/q attraversano la banca band-stop 1.3/1.6/1.9/2.1 Hz prima dei body flap.
- Integrazione e filtri: `dt=0.005 s`, `fs=200 Hz`, Nyquist 100 Hz; filtri e attuatore Tustin sono ricalcolati allo stesso passo.
- Aerodinamica: 58 elementi C81 quasi-stazionari. Otto patch passive delle winglet riproducono area, rastremazione e quarto di corda dei CAERO1 NASTRAN 146001–149001 e 246001–249001.
