# Rapporto: `_trim_coupled_stage`

## Ruolo e stato

È uno staging autonomo per un unico trim accoppiato DUST–MBDyn a 55 m/s, senza
SAS, manovra o BFF. Contiene 31 file e circa 55 MB, quasi tutti dovuti alla
copia `mbdyn_modal_60.fem`. È **preparato ma non chiuso da un risultato locale**:
`results/README.txt` è sostanzialmente vuoto e il README rimanda l'esecuzione al
server.

## Contenuto e codice

- `case/`: FEM a 60 modi, input MBDyn/DUST, mesh, 59 nodi preCICE, XML e patch
  per DUST con preCICE 3/hinge accoppiate;
- `run_trim.py`: risolve i binari, renderizza input/mesh, valida geometria,
  implementa il socket MBDyn, coordina il coupling e valuta residui/pendenze;
- `make_paraview.py`: ripara il VTU raw-appended e crea il PVD;
- `server_instructions.txt` e README: build patchata ed esecuzione remota;
- `results/`: placeholder.

Il caso trattiene heave e pitch rigidi, lascia liberi i modi elastici 7–18 e
muove simmetricamente solo BFL/BFR. I carichi sono rampati, il trim è integrato
online, poi i comandi vengono congelati per la verifica finale.

## Obsolescenza e azioni

Non è obsoleto come esperimento isolato, ma il nome `_stage` segnala un
incubatore. Parte della sua architettura è stata assorbita da `BFF_DUST_55`.
Eseguire/archiviare un risultato valido oppure integrare il caso come test di
trim del workflow coupled e spostare la patch DUST in una posizione canonica.
Solo dopo questa decisione lo staging può passare in archivio.

