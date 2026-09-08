# Verifiche del 2026-09-08

Completati sei test automatici (`python3 -m unittest discover -s tests -v`):
selezione modi, rimozione solo yaw, filtri stabili a guadagno DC unitario,
punti cerniera distinti, conversione velocita' e soglie di rilascio.

`runs/check_20260908_130829/check.json`: preprocessor DUST, validatore preCICE
e parsing MBDyn con FEM25 completati; geometria `[7320,7503,59,10]`.

`runs/smoke_20260908_131340/analysis.json`: test tecnico accoppiato completo
fino a 0.04 s, 21 stati MBDyn, 25 modi, valori finiti, massimo quattro
iterazioni implicite. Uscite solutori `[0,0]`. Violazione massima del vincolo
yaw 7.75e-11 rad. Variazione delle otto superfici congelate non interessate
dal rap inferiore a 6.4e-13 rad. Nove frame DUST e serie VTK generate.

La riga terminale del socket `got ABORT from peer` compare alla chiusura:
non e' stata usata come criterio di successo; sono stati controllati sia i
codici dei due solutori sia il tempo finale e i risultati.

Nel sorgente DUST `src/geo/mod_hinges.f90` il campo HDF5 `theta` e' lasciato
a zero per `input_type=coupling`: la deformazione accoppiata viene applicata
da `src/precice/mod_precice.f90:update_elems`. Non interpretare quel campo
come misura della deflessione dei flap; usare giunti MBDyn e geometria DUST.

Non eseguita la simulazione lunga. Non validati assestamento, recupero fisico,
taratura PID, convergenza temporale o velocita' di flutter DUST. Rimane da
verificare la fase temporale dei carichi del binario/adattatore ereditati,
come descritto nel README. I risultati di identificazione sono candidati,
non una certificazione del confine di flutter.
