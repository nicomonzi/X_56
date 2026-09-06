# Rapporto: `bff_longitudinal`

## Ruolo e stato

Deriva da `bff_eigenalaisy` e tenta di rendere il SAS trasparente nella banda
BFF con quattro notch longitudinali. Contiene 18 file e circa 56 MB, dominati
dal FEM a 60 modi. Il README dichiara esplicitamente che il progetto **non è
stato simulato durante la preparazione**. È quindi un ramo preparatorio
superseded da `BFF_open_loop`.

## Codice e modifiche

- correzione dei bias scambiati WF1L/BFR;
- notch a 1,3/1,6/1,9/2,1 Hz su pitch e q;
- 25 modi elastici, FEM 7–31, `dt=0,005 s`;
- `run_sweep.py` per 57,5–65 m/s;
- `analyse_bff.py` estende l'analizzatore storico con filtri di secondo ordine e
  band-stop longitudinale;
- `INCLUDE/` contiene il modello autonomo e la C81.

## Obsolescenza e azioni

L'idea dei notch è importante nella genealogia, ma non esistono risultati che
la validino. `BFF_open_loop` adotta un protocollo più forte: blocco completo
delle superfici e feedback nullo. Archiviare come design non eseguito,
conservando README e differenze rispetto al predecessore; non presentarlo come
risultato.

