# Rapporto: `bff_eigenalaisy`

## Ruolo e stato

Studio MBDyn con quattro gradi rigidi liberi e aerodinamica strip-C81
closed-loop. Il nome contiene un refuso storico. Ha 28 file e circa 64 MB:
FEM/INCLUDE, una run da 25,4 MB, CSV e figure. È **superseded**.

## Codice e risultati

`run_sweep.py` varia V_INF; `analyse_bff.py` applica Matrix Pencil,
ricostruisce forze e comandi e genera il report. La copia `analyse_bff copy.py`
è identica alle copie presenti in altri due rami. `MODEL_AUDIT.md` documenta
modi 7–12, vincoli X/Y, controlli, notch e assenza di propulsione.

Il solo punto a 60,8421 m/s ha trim valido e BFF stimato a 1,583 Hz con
`sigma=-0,280 1/s`, ma il verdetto globale è “inconclusive” perché il modo
longitudinale rigido non è identificato con sufficiente evidenza. Il SAS resta
attivo e può partecipare al ramo osservato.

## Obsolescenza e azioni

È stato sostituito prima da `bff_longitudinal` e poi da `BFF_open_loop`, che
impone una vera finestra di feedback nullo. Conservare report, CSV e singolo
NetCDF come evidenza; archiviare il resto dopo aver verificato che FEM e include
siano coperti dall'asset canonico.

