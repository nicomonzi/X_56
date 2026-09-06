# Completamento dello sweep di trim fino a 70 m/s

## Avvio

Il comando va eseguito in WSL, perché l'eseguibile MBDyn configurato è Linux.
Il percorso Windows `C:\Users\Utente\Desktop\TRIM` corrisponde in WSL a
`/mnt/c/Users/Utente/Desktop/TRIM`.

```bash
cd /home/nicomonzi/X_56/workflows/trim
python3 run_trim_sweep.py --continue-on-error
```

Le impostazioni predefinite sono già quelle richieste:

- input `main_trim.mbd`;
- risultati in `/mnt/c/Users/Utente/Desktop/TRIM`;
- velocità da 30 a 70 m/s, passo 2.5 m/s;
- densità `9.7284e-8 lbf s^2/in^4`;
- conservazione dei casi completi e lancio dei soli casi mancanti;
- analisi e grafici automatici alla fine.

La cartella dei risultati viene predisposta vuota: il comando esegue quindi
tutti i 17 casi da 30 a 70 m/s. `--continue-on-error` fa proseguire gli altri
casi se uno fallisce. Senza questa opzione lo sweep si ferma al primo errore.
`--force` riesegue e sovrascrive anche eventuali casi completi e va quindi usato
soltanto intenzionalmente.

Per controllare la configurazione senza lanciare MBDyn:

```bash
python3 run_trim_sweep.py --dry-run --no-analyse --output-dir /tmp/trim_check
```

## Cosa fa `run_trim_sweep.py`

1. Trova MBDyn da `--mbdyn`, dalla variabile `MBDYN_BIN`, dal `PATH` oppure nei
   percorsi di installazione noti.
2. Costruisce la griglia 30, 32.5, ..., 70 m/s senza errori di accumulo del
   passo floating-point.
3. Per ogni velocità crea `V_XXX.X_mps` e verifica `result.nc` e il tempo finale
   riportato in `result.out`. Se il caso ha raggiunto 24 s, lo considera
   completo e non modifica neppure il suo `case_input.mbd`.
4. Per un caso mancante copia il template e sostituisce unicamente `V_INF` e
   `RHO_AIR`. La copia viene conservata come prova esatta dell'input usato.
5. Esegue MBDyn dalla cartella `TRIM`, condizione necessaria per risolvere gli
   include relativi. Salva gli output con prefisso `result` e tutta la console
   in `console.log`.
6. Aggiorna `sweep_manifest.csv` dopo ogni caso, così rimane una traccia anche
   in caso di interruzione.
7. A sweep concluso avvia `analyze_trim_sweep.py` (disattivabile con
   `--no-analyse`).

L'analizzatore legge i NetCDF, converte angoli e reazioni in unità SI, salva un
`history.csv` per caso e produce in `analysis`:

- `trim_summary.csv`, con soluzione, residui, pendenze finali, picchi e stato;
- `analysis_report.txt`;
- `01_trim_solution.png`;
- `02_trim_quality.png`;
- `03_selected_histories.png`.

I grafici seguono lo stile tesi: font serif Computer Modern, testi inglesi,
etichette espresse esclusivamente come simbolo fisico e unità, titoli brevi e
non in grassetto sui singoli pannelli, palette scura, tick interni e griglia
principale leggera. Ogni figura viene salvata sia come PNG a 300 dpi sia come
PDF vettoriale. Le storie selezionate sono quelle più vicine a 30, 50 e 70 m/s.

## Che controllo viene usato per il trim

Non è un PID di volo classico. È un risolutore dinamico del trim longitudinale
con due residui e due incognite:

```text
r = [Fz, My]^T
u = [pitch, elevator]^T
```

Vicino al trim, `r` può essere linearizzato come `r ~= J (u-u_trim)`, dove `J`
è il Jacobiano 2x2 delle reazioni rispetto ai due comandi. Le combinazioni
`TRIM_JINV_ij` applicano una stima di `J^-1` ai residui, disaccoppiando i due
canali. Ciascun blocco PID usa poi `Kp=0`, `Ki=TRIM_RESPONSE_RATE=0.12` e
`Kd=0`: è quindi un puro integratore. L'integrale è adatto al problema statico
perché continua a correggere finché forza verticale e momento di beccheggio
sono nulli. I valori `Ii0` sono stime iniziali di pitch ed elevatore: riducono il
transitorio, ma non dovrebbero determinare la soluzione finale.

Il Jacobiano è stato identificato alle condizioni di riferimento e il suo
inverso viene scalato con

```text
(rho_ref/rho) (V_ref/V)^2
```

perché, in prima approssimazione, le derivate aerodinamiche dimensionali sono
proporzionali alla pressione dinamica `q = 0.5 rho V^2`.

## Come ricavare e tarare i guadagni

Per ricavare la matrice di disaccoppiamento a una velocità di riferimento:

1. scegliere un punto vicino al trim e disattivare temporaneamente
   l'integratore;
2. perturbare il pitch di `+dtheta` e `-dtheta`, lasciando fisso l'elevatore;
3. perturbare l'elevatore di `+ddelta` e `-ddelta`, lasciando fisso il pitch;
4. dopo il transitorio, misurare le medie di `Fz` e `My`;
5. calcolare con differenze centrate
   `dFz/dtheta`, `dMy/dtheta`, `dFz/ddelta`, `dMy/ddelta`;
6. assemblare `J`, controllarne numero di condizionamento e segni, quindi
   calcolare `J^-1` e inserirlo nei quattro `TRIM_JINV_ij`;
7. validare con perturbazioni diverse da quelle usate nell'identificazione.

Con un disaccoppiamento ideale, la dinamica locale del puro integratore è
approssimativamente esponenziale e `Ki` si comporta come una velocità di
convergenza in `1/s`: la costante di tempo è circa `tau = 1/Ki`. Con
`Ki=0.12`, `tau` è circa 8.3 s. Il controllo è attivo approssimativamente da 2
a 20 s; senza l'aiuto del seed, in 18 s una dinamica ideale con questo guadagno
eliminerebbe circa l'88% dell'errore iniziale.

Una taratura pratica consiste nel provare, per esempio, `Ki = 0.08, 0.12,
0.16, 0.20, 0.24`, mantenendo invariati filtro e Jacobiano, e scegliere il
valore più alto che:

- non genera oscillazioni lente o eccitazione dei modi elastici;
- non porta pitch/elevatore ai limiti di saturazione;
- rende piccoli `Fz_mean_N` e `My_mean_Nm` dopo il congelamento;
- rende quasi nulle le pendenze dei comandi fra 18 e 20 s;
- resta robusto su tutta la griglia di velocità.

Se si desidera almeno il 95% di riduzione ideale entro 18 s, un punto di
partenza teorico è `Ki ~= 3/18 = 0.167 1/s`; per il 99% è circa 0.256 1/s. Sono
solo punti di partenza: filtro, saturazioni, non linearità e dinamica elastica
impongono la verifica numerica. Per questo problema MIMO, una taratura tipo
Ziegler--Nichols è meno indicata della stima del Jacobiano e di uno sweep del
guadagno.

Il Jacobiano è stato re-identificato a 63 m/s e densità `1.146e-7` IPS usando
esclusivamente BFL/BFR. Sono state applicate differenze finite centrate di
0.25 deg sul pitch e 0.50 deg sui body flap, con `Ki=0`; le reazioni sono state
mediate fra 22 e 24 s. Input, NetCDF, tabella e report riproducibile sono in
`jacobian_calibration_bfl_bfr`. La validazione indipendente del passo Newton ha
ridotto i residui da `[-71.10 lbf, 631.29 lbf in]` a
`[1.52 lbf, 1.10 lbf in]`, inferiori alla deviazione standard finale. I nuovi
coefficienti `TRIM_JINV_ij` e i seed alla condizione di riferimento sono già
inseriti in `main_trim.mbd`.

## Perché serve il filtro Butterworth

`Fz` e `My` sono reazioni istantanee di un modello flessibile. Contengono il
valore quasi statico necessario al trim, ma anche oscillazioni modali, rumore
numerico e il transitorio iniziale. Se questi contributi entrano direttamente
nell'integratore, il comando può inseguire le vibrazioni, accumulare errore
inutile, oscillare o raggiungere la saturazione.

Il filtro è un passa-basso Butterworth digitale del secondo ordine, con
`fc=0.25 Hz`, `dt=0.02 s` e frequenza di campionamento `fs=50 Hz`. Butterworth è
adatto perché ha banda passante massimamente piatta: non introduce ripple nel
segnale quasi statico. La sua equazione è

```text
y[k] = A1 y[k-1] + A2 y[k-2]
     + B0 x[k] + B1 x[k-1] + B2 x[k-2]
```

con i coefficienti `BW_A1`, `BW_A2`, `BW_B0`, `BW_B1`, `BW_B2` presenti nel
file MBDyn. Per ricalcolarli in Python:

```python
from scipy.signal import butter

dt = 0.02
fc = 0.25
fs = 1.0 / dt
b, a = butter(2, fc / (fs / 2.0), btype="low")
mbdyn_a1, mbdyn_a2 = -a[1], -a[2]
mbdyn_b0, mbdyn_b1, mbdyn_b2 = b
```

Una frequenza di taglio più alta rende il trim più reattivo ma lascia passare
più contenuto elastico; una più bassa pulisce meglio il segnale ma aggiunge
ritardo e può impedire la convergenza prima del freeze. La rampa `step5` fra 1
e 2 s completa la protezione evitando che il controllore riceva bruscamente
l'intera reazione iniziale.
