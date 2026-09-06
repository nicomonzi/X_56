# X-56 BFF durante dive--pull-up

Pipeline MBDyn definitiva per misurare come una manovra longitudinale modifica
l'onset del body-freedom flutter nell'intorno della frontiera steady di circa
66.03 m/s.

## Campagna

- velocita': 65, 66 e 67 m/s;
- classi nominali di aggressivita': 1.15, 1.30, 1.60 e 1.80;
- una coppia `shadow/excited` per punto;
- 24 traiettorie complessive;
- passo temporale 0.02 s;
- SAS-off lungo 2.05 s.

Le classi nominali non sono setpoint del fattore di carico. Servono soltanto a
ricavare il comando cinematico di pitch-rate. Non esiste feedback su `n`: il
carico medio, minimo, massimo, deviazione e deriva vengono ricostruiti a
posteriori dalla shadow.

La manovra parte livellata, entra in una picchiata controllata, compie il
pull-up e torna al trim. Profondita' della picchiata e durata del recupero sono
adattate automaticamente all'aggressivita'. Durante il SAS-off tutte le
superfici sono sample-and-hold al valore di rilascio; soltanto la run excited
aggiunge il piccolo rap simmetrico sulle WF4.

## Esecuzione

Controllare il piano senza avviare MBDyn:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 run_dive_pullup_sweep.py
```

Prima fase consigliata, quattro shadow a 66 m/s:

```bash
python3 run_dive_pullup_sweep.py \
  --velocities 66 \
  --shadow-only \
  --execute \
  --jobs 1
```

Campagna completa:

```bash
python3 run_dive_pullup_sweep.py --execute --jobs 1
```

Le shadow definitive gia' presenti vengono riutilizzate automaticamente. I
risultati sono salvati in:

```text
C:\Users\Utente\Desktop\BBF_PULLUP\dive_pullup_sweep
```

## Analisi

Al termine della campagna lo script esegue automaticamente:

```bash
python3 analyse_dive_pullup.py \
  /mnt/c/Users/Utente/Desktop/BBF_PULLUP/dive_pullup_sweep
```

L'analisi usa `excited-shadow`, verifica modo 7 e deformazione simmetrica della
tip, controlla il leakage delle superfici e genera CSV, JSON e grafici. Un
punto viene escluso se il carico non e' quasi stazionario durante il SAS-off,
modo e tip non concordano o la traiettoria non torna quasi livellata.

## Limite del modello

Il modal joint conserva massa, smorzamento e rigidezza del ROM corrente. I
risultati misurano gli effetti cinematici, aerodinamici e dei controlli della
manovra; non includono lo stress stiffening dipendente dal precarico.
