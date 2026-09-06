# Audit del modello

- Nodo modale: `990000`; giunto modale: `5`. Il total-pin `1` ha due componenti attive: traslazioni globali X e Y.
- Gradi rigidi liberi: Z, roll, pitch e yaw. X e Y sono le sole componenti traslazionali `active`; Z e tutte le rotazioni del joint `1` sono `inactive`.
- Modi FEM attivi: 7–12. Modi 1–6: rigidi esclusi. Frequenze secche: 3.2171, 5.3027, 8.7051, 11.1640, 12.2571, 12.7589 Hz.
- Superfici: BFL/BFR `1004/2004`; WF1 `1008/2008`; WF2 `1011/2011`; WF3 `1014/2014`; WF4 `1017/2017`.
- Combinazione simmetrica: stesso segno numerico L/R. Combinazione differenziale: segno opposto L/R.
- Rappresentazione del moto: velocità strutturale iniziale nulla e vento uniforme `+VINF`; `Vrel` viene sempre ricalcolata come vento meno velocità inerziale.
- Controllo: quota/Vz → riferimento pitch; pitch/q → body flap simmetrici; roll/p → WF1/WF2 differenziali; yaw/r/Vy → body flap differenziali.
- Non esistono propulsione, throttle, airspeed hold o feedback delle reazioni X/Y.
- Il trim e' un termine separato; le correzioni degli attuatori partono da zero e non esiste alcun rilascio programmato del SAS.
- Burst: WF4 simmetrico, finestra Hann, 4 cicli a 2.0597 Hz. Il SAS resta attivo e notchato durante e dopo il burst.
- Aerodinamica: 58 elementi C81 quasi-stazionari. Otto patch passive delle winglet riproducono area, rastremazione e quarto di corda dei CAERO1 NASTRAN 146001–149001 e 246001–249001.
