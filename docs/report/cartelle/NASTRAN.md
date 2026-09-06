# Rapporto: `NASTRAN`

## Ruolo e stato

È il deposito principale dei modelli Nastran originali e rilasciati: 8.885 file
e circa 559 MB. Il conteggio è gonfiato da un ambiente Python completo
versionato dentro `REALASED_MODEL/MAIN/.venv_x56`: 3.530 `.py`, 3.527 `.pyc`,
154 `.so` e centinaia di stub. La cartella è **riferimento misto**, non tutta
obsoleta.

## Sottocartelle

- `NASTRAN40` e `FEMGEN40`: sviluppo della base modale e conversione FEM;
- `NASTRAN_ALTER_OK`: variante con ALTER verificata;
- `SOL101`: modello statico e output;
- `FEMParser_OK`: parser del formato FEM MBDyn e file di prova;
- `REALASED_MODEL` (nome con refuso conservato): modello rilasciato, con `BULK`,
  `FLUTTER_TEST`, `TRIM`, `Airfoil` e `MAIN`. Gli script visualizzano airfoil,
  forme modali e generano VTU/PVD.

## Codice

`FEMParser.py` legge header, nodi e matrici del FEM. Dentro il modello rilasciato
`plot_airfoil_sections.py` produce tavole spanwise; `plot_x56_modes.py` e
`view_x56_modes.py` usano pyNastran per plot statici/interattivi;
`export_x56_paraview.py` esporta geometria e modi in VTU/PVD. Gli shell script
installano/avviano il viewer.

## Problemi di struttura

`.venv_x56` occupa circa 348 MB ed è obsoleto/rigenerabile; da solo spiega la
maggior parte dei file. Sono presenti 3.527 bytecode `.pyc` dentro il venv e
molti pacchetti terzi che non dovrebbero essere versionati. I deck BULK sono
duplicati in almeno otto workflow, ma rappresentano una baseline condivisa.

## Azioni

Mantenere `REALASED_MODEL`, i deck, F06/OP2/MAT di riferimento e i tool scritti
per il progetto. Rimuovere dal versionamento l'ambiente virtuale solo con
commit dedicato, sostituendolo con `requirements.txt`/lock e istruzioni di
creazione. Centralizzare BULK per versione del modello; lasciare nei casi un
manifest o uno script di materializzazione. Rinominare `REALASED_MODEL` solo in
una migrazione controllata perché molti path possono dipenderne.

