# Rapporto: `GRAFICI`

## Ruolo e stato

Raccoglie figure e workflow di confronto destinati alla tesi. Ha 349 file e
circa 150 MB. È un **deliverable attivo**, ma mescola sorgenti, copie dei modelli,
output solver e figure finali.

## Sottocartelle e codice

- `MODAL_COMPARISON_10LB/`: `compare_mode_shapes.py` estrae forme da FEM/F06 e
  MOV, allinea segno/normalizzazione, calcola MAC e confronta i modi;
  `generate_mbdyn_comparisons.py` esegue i modi 7–18;
  `interactive_modal_plot.py` offre un viewer 3-D. Include copie autonome
  Nastran/MBDyn e risultati.
- `SOL101/`: `compare_deformed_shapes.py` confronta deformate statiche per tre
  load case, con modelli Nastran/MBDyn locali.
- `X56_BEAM_GEOMETRY_PLOTS/`: parser BDF e figure 2-D/3-D della geometria
  strutturale e della trave fittizia.
- `X56_PARAVIEW/`: 60 MB di VTU/PVD e immagini per la visualizzazione modale.

## Obsolescenza e azioni

Non è obsoleta come prodotto, ma le copie BULK/FEM e gli output solver non
dovrebbero vivere insieme alle tavole finali. Separare `scripts/`, `data/` e
`figures/`; registrare per ogni figura comando, input hash e commit. I 15 file
extra duplicati interni occupano poco, ma molti duplicati attraversano
`NASTRAN` e `NASTRAN_SIMULATIONS`: mantenerli solo se servono alla riproducibilità
standalone.

