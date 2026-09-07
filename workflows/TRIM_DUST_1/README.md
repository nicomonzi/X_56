# X-56A: trim DUST–MBDyn a 55 m/s

Il caso esegue **un solo accoppiamento** DUST–MBDyn. Non contiene SAS,
manovre, rilascio del velivolo o sweep.

## Modello

- Mesh validata di `DUST_MESH`: 2440 pannelli, 2583 punti, nessuna winglet.
- Passo fisso `dt = 0.004 s`; modi elastici FEM 7--18, fino a 18.559 Hz.
- FEM copiato: `case/mbdyn_modal_60.fem` (60 modi disponibili).
- 59 nodi preCICE: 39 nodi elastici e 20 estremi delle dieci cerniere.
- Tutte le hinge usano `hinge_rotation_input = coupling`. DUST riceve quindi
  posizione, rotazione, velocità e velocità angolare dei veri assi MBDyn.
- Nel trim si muovono simmetricamente solo BFL e BFR; WF1--WF4 sono a zero ma
  restano già collegati per il successivo SAS/BFF.
- Heave e pitch rigidi sono trattenuti. I modi elastici restano liberi.

## Scaletta temporale

- `0.0--0.4 s`: rampa dolce dei carichi DUST;
- `0.4--0.8 s`: inserimento dolce dei due integratori di trim;
- `0.4--6.0 s`: pitch e BFL/BFR vengono corretti online nella stessa run;
- `6.0--8.0 s`: i due comandi sono congelati e si verifica l'assestamento.

Il risultato è accettato soltanto se, nella finestra congelata, le reazioni
del giunto di trim e le loro pendenze sono entro le tolleranze. Incidenza,
BFL/BFR, residui e pendenze vengono salvati in `results/trim_55.json`.

## Patch pulita della build DUST sul server

La patch parte dal commit DUST originale del server e comprende sia la
chiusura corretta con preCICE 3 sia l'aggiornamento della griglia virtuale per
le hinge accoppiate. Prima si ripristinano **solo i due sorgenti interessati**:

```bash
cd ~/dust-patched

cp /home/dust-group/software/dust/src/dust.f90 src/dust.f90
cp /home/dust-group/software/dust/src/precice/mod_precice.f90 \
   src/precice/mod_precice.f90

patch --dry-run -p1 \
  < ~/ZENO/TRIM_55_DUST/case/dust_precice3_coupled_hinges.patch
patch -p1 \
  < ~/ZENO/TRIM_55_DUST/case/dust_precice3_coupled_hinges.patch

cmake --build build-user -j24
```

Non usare la build standalone `~/.local/dust-0.8.2-b`: non contiene preCICE.

## Esecuzione

```bash
cd ~/ZENO/TRIM_55_DUST
unset DISPLAY
export DUST_BIN="$HOME/dust-patched/build-user/bin/dust"
export DUST_PRE_BIN="$HOME/dust-patched/build-user/bin/dust_pre"

python3 run_trim.py --check --threads 24
python3 run_trim.py --threads 24
```

Monitoraggio:

```bash
tail -F ~/ZENO/TRIM_55_DUST/work/current/coupling.csv
```

Non aggiungere `work/` a Git: contiene gli HDF5 e gli output della run.

## ParaView

Al termine:

```bash
cd ~/ZENO/TRIM_55_DUST/work/current
mkdir -p paraview
"$HOME/dust-patched/build-user/bin/dust_post" ../../case/dust_post.in
cd ../..
python3 make_paraview.py
paraview work/current/paraview/trim55_all.pvd
```

Il PVD usa le coordinate scritte da DUST: deformazione elastica, superfici
mobili e scia appartengono allo stesso riferimento e alla stessa timeline.
