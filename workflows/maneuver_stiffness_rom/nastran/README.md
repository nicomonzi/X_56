# Parte Nastran: SOL 103 con precarico statico

`prepare_cases.py` prepara due casi autosufficienti a 1.0 g e 1.6 g partendo
dai deck SOL 103/`STATSUB=1` che hanno già completato correttamente. Non avvia
Nastran. I file eseguibili vengono scritti in
`/home/nicomonzi/ZENO/prestress_stiffness_rom`.

```bash
python3 nastran/prepare_cases.py --overwrite
```

L'opzione `--overwrite` elimina la vecchia campagna `linear_qrmeth3`, i deck e
gli output SOL 106 falliti e gli eventuali output omonimi SOL 103. Conserva e
rigenera carichi, RBE3 e cartella `BULK` dalle sorgenti validate.

## Run manuali su Zeno

Sono richiesti soltanto questi due job:

```bash
cd /home/nicomonzi/ZENO/prestress_stiffness_rom/n1p000
./RUN_BY_USER.sh

cd /home/nicomonzi/ZENO/prestress_stiffness_rom/n1p600
./RUN_BY_USER.sh
```

Se il comando del solver non è `nast`:

```bash
NASTRAN_CMD=nastran ./RUN_BY_USER.sh
```

Ogni job deve produrre tre file non vuoti:

- `prestress_sol103_statsub.f06`;
- `prestress_sol103_statsub.op2`;
- `mbdyn_modal.mat`.

Il runner rimuove i propri output prima del lancio e controlla fatal error,
`END OF JOB` e presenza dei tre risultati.

## Formulazione

Il subcase 1 calcola l'equilibrio statico lineare con `LOAD=613` e `SPC=620`.
Il subcase 2 usa `STATSUB=1` per includere la rigidezza differenziale, ma non
ha un SPC: i modi restano free-free come la base MBDyn. `PARAM,FOLLOWK,NO`
rende esplicita l'esclusione della rigidezza da forze follower, coerente con
lo snapshot dei carichi nel riferimento body.

Si mantiene la formulazione originale CQUADR/CTRIAR, senza `QRMETH=3`. L'ALTER
`MBDyn_NASTRAN_alter_SOL103_2024.nas` esporta `MHH`, `KHH` e `LUMPMS`; OP2
contiene geometria e forme modali su tutti i GRID. Questo è necessario per
proiettare sulla `Phi0` MBDyn completa, che contiene 8527 nodi.

`KHH` è espresso nella base modale del singolo snapshot. Dopo i run si usa
`femgen`, quindi `project_fixed_basis.py`, per ricostruire

```text
Phi0 ~= Psi_n T_n
K_h^Phi0(n) = T_n^T KHH_n T_n
```

nella base fissa dei modi FEM 7--12. Il programma misura il residuo della
ricostruzione e rifiuta una proiezione oltre la soglia configurata.
