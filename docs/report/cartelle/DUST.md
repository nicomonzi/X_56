# Rapporto: `DUST`

## Ruolo e contenuto

Piccola cartella sorgente, 15 file e circa 21 kB. Contiene gli input
`parametric_mesh.in`, `dust_pre.in`, `dust.in`, `dustPost_IntLoads.in` e nove
sezioni aerodinamiche `.dat` sotto `airfoilsection/`. `convert_paraview.py`
legge output DUST raw-appended, li converte e costruisce una collezione PVD.
`Output/` e `Postpro/` sono vuote nello snapshot.

La geometria parametrica e gli airfoil sono l'origine dichiarata del primo
smoke in `BFF_DUST_55`; il lavoro di luglio ha aggiunto le superfici mobili.
Non contiene winglet e questo limite si propaga ai casi accoppiati recenti.

## Stato e obsolescenza

È **supporto corrente**, non obsoleta. Il problema è la duplicazione: copie
delle sezioni compaiono in molti frame di animazione e casi DUST. Conviene
promuovere qui, o in `models/dust/x56`, una versione canonica con SHA-256 e
fare generare le copie runtime. Non cancellare finché i runner recenti non sono
stati aggiornati e verificati.

