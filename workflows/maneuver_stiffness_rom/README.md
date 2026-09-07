# X-56 maneuver stiffness ROM

L'integrazione operativa nella manovra BFF è la campagna
`../maneuver_bff/run_sweep.py --campaign prestress_rom`. Usa la manovra
`dive_pullup` validata, conserva la finestra SAS-off e programma
`Delta n=VINF*q_command/GRAVITY`; non usare il vecchio pull-up con SAS continuo.

Questo workflow implementa la strada incrementale controllabile:

```text
K_h(n_z) = K_h(1g) + (n_z - 1) K_h,n
Q_prestress = -(n_z - 1) K_h,n (q - q_eq)
```

## Stato dell'implementazione

I due snapshot SOL103/STATSUB free-free a 1.0 g e 1.6 g sono stati convertiti
con `femgen` e proiettati nella base MBDyn. I residui massimi di sottospazio
sono rispettivamente 0.088% e 0.142%; il mismatch tra la matrice proiettata a
1 g e il `.fem` baseline è 0.152%. La ROM e la forza modale MBDyn sono in
`generated/` e un caso 1.6 g ha superato parsing e avanzamento temporale.

Con due soli punti la secante 1.0--1.6 g è definita, ma la linearità rispetto a
`n_z` non è ancora verificata indipendentemente; per farlo servirà in seguito
uno snapshot intermedio, per esempio a 1.3 g.

La base modale resta quella del `.fem` validato, con i modi FEM 7--12. Non si
modificano `femgen`, MBDyn o il `.fem`; il termine aggiuntivo e una forza modale
matriciale 6 x 6.

## Precisazione fisica fondamentale

Il `.fem` attuale e il modello dry/free-free baseline, non una matrice ottenuta
da un vero equilibrio Nastran prestressato a 1 g. Di conseguenza la prima
implementazione e rigorosamente una **correzione incrementale che si annulla a
1 g**:

```text
K_eff(n_z) = K_fem + [K_h(n_z) - K_h(1g)].
```

Essa coincide con `K_eff = K_h(n_z)` soltanto se `K_fem = K_h(1g)`. Il report
generato misura questa differenza e la segnala; non la nasconde.

## Procedura

1. I run sono preparati in `ZENO`. Per rigenerarli, senza eseguirli:

   ```bash
   cd /home/nicomonzi/X_56/workflows/maneuver_stiffness_rom
   python3 nastran/prepare_cases.py --overwrite
   ```

2. Eseguire manualmente i due casi SOL 103 con `STATSUB=1` nelle cartelle
   `/home/nicomonzi/ZENO/prestress_stiffness_rom/n1p000` e `n1p600`. Entrambi
   mantengono gli elementi CQUADR/CTRIAR originali, senza `QRMETH=3`, e usano
   `FOLLOWK=NO`. Il precarico è supportato nel subcase statico, mentre il
   subcase modale resta free-free come il `.fem` MBDyn. I dettagli sono in
   `nastran/README.md`.

3. Dopo avermi restituito gli output, convertire ogni snapshot con `femgen` e
   proiettarlo nella base fissa. Per esempio per 1.0 g:

   ```bash
   cd /home/nicomonzi/ZENO/prestress_stiffness_rom/n1p000
   /usr/local/mbdyn/bin/femgen prestress_sol103_statsub \
       -o prestress_sol103_statsub.fem
   cd /home/nicomonzi/X_56/workflows/maneuver_stiffness_rom
   python3 project_fixed_basis.py \
       --case-fem /home/nicomonzi/ZENO/prestress_stiffness_rom/n1p000/prestress_sol103_statsub.fem \
       --load-factor 1.0 --output matrices/n1p0
   ```

   Si ripete con `n1p600`, `--load-factor 1.6` e `matrices/n1p6`. La proiezione
   usa tutte le forme prestressate, misura il residuo di sottospazio e rifiuta
   il risultato oltre il 5%.

4. Costruire la ROM:

   ```bash
   python3 prestress_rom.py
   ```

   Il fit e vincolato a passare per il caso 1 g. Con soli due punti produce la
   secante, ma marca `linearity_validated=false`; un terzo punto (per esempio
   1.3 g) e necessario per verificare davvero l'ipotesi lineare.

5. Generare direttamente il caso di manovra con la ROM attiva:

   ```bash
   python3 ../maneuver_bff/maneuver_case.py --family pullup --velocity 66.75 \
       --load-factor 1.6 --shadow --prestress-rom --sas-continuous --dry-run \
       --output results/prestress
   ```

   Per eseguirlo basta rimuovere `--dry-run`. Il generatore conserva la
   correzione DLM esistente, aumenta il numero di forze e aggiunge
   `PRESTRESS_ROM_FORCE`. Il suffisso `_prestress_commanded_sas_continuous`
   impedisce di sovrascrivere il caso baseline. Senza `--sas-continuous` il
   caso torna al protocollo di identificazione BFF open-loop e può divergere
   fisicamente. `render_mbdyn_case.py` resta disponibile
   per convertire manualmente casi `.mbd` già generati.

## Scelta di n_z

Il default `commanded` usa la classe di carico della manovra e il relativo
profilo temporale. E la scelta raccomandata per il primo confronto A/B perche
replica esattamente la parametrizzazione dei casi Nastran e non introduce un
feedback da accelerazione. Il renderer la accetta solo per la famiglia
`pullup`; per `pullup_angle`, `dive_pullup` e `roll` occorre una legge dedicata
oppure `measured`.

L'opzione

```bash
python3 render_mbdyn_case.py ... --nz-source measured
```

usa invece il carico specifico body-z istantaneo:

```text
n_z = xPP[3]/g + cos(roll) cos(pitch).
```

Va trattata come analisi di sensibilita: la forza modale standard basata su
drive e esplicita e non aggiunge un Jacobiano strutturale; usare `xPP` crea
quindi un anello accelerazione-forza che puo peggiorare la convergenza.

## q_eq e confronti minimi

La prima versione usa i `Q_INIT_7...Q_INIT_12` del trim schedulato. Le prove da
eseguire in coppia, con stesso stato e stessa eccitazione, sono:

- ROM disattiva contro ROM attiva a 1 g: differenza numerica nulla;
- ROM disattiva contro ROM attiva a 1.6 g;
- shadow contro excited per non attribuire alla manovra il transitorio BFF;
- `dt=0.01` contro `dt=0.005` al punto critico;
- controllo che gli autovalori di `K_eff` restino fisicamente accettabili.

Il parametro `n_z` resta una riduzione dello stato di tensione, non una sua
descrizione univoca. La validita ottenuta vale per la famiglia di manovre e la
distribuzione di carico usate nei run Nastran. In questa prima campagna la
coordinata del fit e la classe nominale 1.0/1.6; nei manifest sono conservati
anche i valori body-z recuperati, 0.99398/1.63290.

## Verifica locale

```bash
python3 -m unittest discover -s tests -v
python3 prestress_rom.py  # fallisce intenzionalmente finche mancano le matrici
```

Riferimenti locali verificati: ALTER SOL103 moderno in
`/home/nicomonzi/src/mbdyn/etc/modal.d/`, esempio shell-blade in
`/home/nicomonzi/src/mbdyn/tests/joints/modal/shell-blade/`, forza modale
esistente in `workflows/bff_open_loop/main_bff_open_loop.mbd`.
