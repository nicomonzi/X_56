# X-56 BFF longitudinale con SAS trasparente

Questo caso deriva da `bff_eigenalaisy`, ma evita che il ramo BFF identificato
sia dominato dal feedback longitudinale del SAS.

## Modifiche principali

- `WF1L` usa `TRIM_SURFACE` e `BFR` usa `BODY_TRIM_SURFACE`: sono corretti i
  due bias scambiati presenti nel caso sorgente.
- I segnali `PITCH_PID` e `Q_PID` attraversano quattro notch in cascata,
  centrati a 1.3, 1.6, 1.9 e 2.1 Hz (`fs=200 Hz`, `Q=2`).
- La cascata attenua almeno circa 20 dB tra 1.25 e 2.2 Hz e circa 36.5 dB a
  2.0597 Hz. Il SAS resta attivo a bassa frequenza per mantenere il moto rigido.
- Il burst simmetrico su WF4 non attraversa il SAS e resta a 2.0597 Hz.
- `INCLUDE/mbdyn_modal.fem` e' la copia locale del FEM SOL103 a 60 modi; il
  modal joint usa i primi 25 modi elastici, cioe' FEM 7--31.
- Il passo e' 0.005 s (`fs=200 Hz`, Nyquist 100 Hz), coerente con il modo 31
  a 36.970 Hz. Tutti i filtri digitali e l'attuatore sono ricalcolati per
  questo passo.
- Nessuna spinta, forza artificiale o vincolo aggiuntivo e' stato introdotto.

Questo progetto non e' stato simulato durante la preparazione. La trasparenza
del SAS deve essere verificata a posteriori confrontando, alla frequenza del
modo identificato, l'ampiezza dei body flap con pitch, q e modo FEM 7.

## Sweep

Da questa cartella:

```bash
python3 run_sweep.py
```

Il comando usa per default `57.5, 60.0, 62.5, 65.0 m/s` e salva i file piatti
in `C:\Users\Utente\Desktop\bbf_longitudinal`, che in WSL corrisponde a
`/mnt/c/Users/Utente/Desktop/bbf_longitudinal`. Per rifare casi esistenti:

```bash
python3 run_sweep.py --overwrite
```

Per scegliere un'altra destinazione:

```bash
BFF_OUTPUT_DIR=/percorso/output python3 run_sweep.py
```

L'analisi resta disponibile con:

```bash
python3 analyse_bff.py
```
