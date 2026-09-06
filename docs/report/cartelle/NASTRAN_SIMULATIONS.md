# Rapporto: `NASTRAN_SIMULATIONS`

## Ruolo e stato

È il pacchetto ordinato di validazione strutturale. Contiene 643 file e circa
455 MB. È **supporto attivo**: la convergenza interna è chiusa, il confronto
esterno Nastran–MBDyn resta aperto.

## Sottocartelle

### `01_SOL103_60_MODES`

Contiene deck, BULK, OP2/F06/MAT e `mbdyn_modal_60.fem`. Lo script
`analyze_modal_spectrum.py` estrae 60 modi: sei quasi rigidi e 54 elastici fra
3,217 e 62,861 Hz. Il FEM ha 8.527 nodi ed è la base da 57,4 MB replicata in
diversi workflow recenti.

### `02_COUPMASS_STUDY`

Confronta massa lumped (`COUPMASS=-1`) e coupled (`COUPMASS=1`).
`compare_coupmass.py` usa MAC sui 45 nodi d'interfaccia e finestra di frequenza;
45 coppie su 54 superano MAC 0,90. Lo scarto medio affidabile è 0,200%, massimo
0,826%. `plot_frequency_comparison.py` produce le figure.

### `03_GRAVITY_5G`

Contiene SOL 101 Nastran a 5 g, template MBDyn e 54 risultati modali.
`run_full_convergence.py` orchestra femgen, verifica FEM, casi e analisi;
`run_remaining_54_bases.py` riprende i mancanti. `analyze_convergence.py`
confronta tip e semispan; script separati calcolano errore locale, cumulativo e
integrale, mentre `plot_modal_truncation_convergence.py` produce il grafico di
tesi.

## Risultati

La prima base sotto 0,01 in rispetto a 54 modi è 25, ma 26 e 36 risalgono sopra
soglia. Da 37 in poi tutte le basi restano sotto; 30 modi è la raccomandazione
operativa robusta adottata dal report. Nessuna base soddisfa però la soglia
esterna: con 54 modi l'errore vettoriale alle tip è circa 1,17 in, 8,4–8,5%.

Il plateau indica un'incoerenza di modello, non insufficienza di modi. Il FEM
SOL 103 è estratto con RBE3 alla radice, il SOL 101 usa RBE2; inoltre una
deformazione di circa 14 in rende delicato il confronto con cinematica lineare.

## Obsolescenza e azioni

La cartella non è obsoleta. Sono rigenerabili `__pycache__`, molti `.log/.out`
e grafici, ma F06/OP2/MAT/FEM e CSV finali sono evidenze. Ogni caso possiede una
copia completa di BULK: utile per autonomia, costosa in duplicati. Prossimo
passo: estrarre un nuovo FEM con lo stesso RBE2 del SOL 101 e ripetere; se il
plateau persiste, usare una soluzione geometricamente non lineare o un carico
più basso.

