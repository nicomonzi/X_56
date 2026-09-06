# X-56 BFF durante manovra

Il modello parte con coordinate modali e superfici nulle. L'unico seed e' una
stima dell'assetto rigido (0.34 deg a 60.8421 m/s), che non costituisce una
deformata elastica. Da `t=0` il sistema porta autonomamente il velivolo
all'equilibrio; la manovra inizia a 14 s.
X e Y restano vincolati, mentre Z, roll, pitch e yaw sono liberi.

## Selezione del caso

```bash
cd /home/nicomonzi/TESI/bbf_manouver

# BFF longitudinale dopo il trim, senza manovra
python3 run_manouver.py --mode level --bff --velocity 60.8421

# Solo pull-up
python3 run_manouver.py --mode pullup --no-bff --velocity 60.8421

# BFF durante il pull-up
python3 run_manouver.py --mode pullup --bff --velocity 60.8421

# Roll con BFF simultaneo
python3 run_manouver.py --mode roll --bff --velocity 60.8421
```

Gli output vengono scritti in `bbf_manouver/output`. Il file base non viene
riscritto: per ogni run viene generato un `.mbd` con nome descrittivo.

## Sequenza temporale

- 0--12 s: acquisizione automatica del trim con vincolo temporaneo;
- 12--14 s: rilascio e verifica del volo libero;
- 14--16 s: ingresso graduale nella manovra;
- 16--20 s: mantenimento;
- 20--22 s: rientro;
- 22--28 s: osservazione del recupero.

Il burst, se abilitato, inizia a 16.25 s e cade nella fase mantenuta della
manovra. Il SAS non viene mai disattivato.

## Controllo

- trim temporaneamente vincolato: integratori su Fz, Mx, My e Mz;
- rilascio a 12 s: restano vincolate soltanto X e Y;
- pull-up: riferimento morbido di pitch -2 deg e Vz +1 m/s;
- roll: riferimento morbido di roll 5 deg;
- burst e manovra sono selezionabili indipendentemente.

I casi validati a 60.8421 m/s sono `level --no-bff`, `pullup --no-bff`,
`roll --no-bff` e `pullup --bff`.
