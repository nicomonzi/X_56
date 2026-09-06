# Rapporto: `X56_AERO_POLAR`

## Ruolo e stato

Studio di confronto aerodinamico Nastran–MBDyn: 131 file e circa 126 MB. Il peso
principale è `nastran/x56_polar.f06` da 70 MB e le copie FEM. È **supporto
storico parzialmente superato**, non un workflow BFF finale.

## Struttura e codice

- `nastran/`: deck polar, BULK e F06;
- `mbdyn/`: modelli e `INCLUDE` con C81 root/wing/winglet e aerobody generati;
- `test/`: `extract_nastran_aero.py` estrae coefficienti dal F06;
  `create_c81_file.py` e `create_spanwise_c81_files.py` costruiscono le polar;
  `create_aerobody_files.py` assegna le sezioni;
  `create_mbdyn_models.py` genera i modelli;
  `run_mbdyn_simulation.py` esegue e raccoglie i risultati;
  `compare_results.py` produce statistiche e grafici.

Il report storico dichiara errore RMS circa 5% e bias medio 0,6%, usando fattori
spanwise 1,15/1,00/0,85. Questi fattori sono una strategia modellistica, non
una distribuzione DLM identificata; vanno mantenuti distinti dalla C81
effettiva e dalla correzione ROM di `BFF_open_loop`.

## Obsolescenza e azioni

Non eliminare il F06 né gli script: sono la provenienza della calibrazione. I
modelli generati e alcune copie `aerobody_orig.mbd`/`aerobody.mbd` sono
rigenerabili. Marcare chiaramente lo studio come confronto preliminare e
spostarlo sotto `validation/aero_polar`; indicare quale C81 è canonica per ogni
workflow per evitare che scaling diversi vengano combinati accidentalmente.

