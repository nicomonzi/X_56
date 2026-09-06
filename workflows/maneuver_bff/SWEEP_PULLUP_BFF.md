# Sweep minimo: il pull-up eccita il BFF?

## Stato e obiettivo

Solo preparazione. Le simulazioni verranno lanciate dall'utente. Non serve
Nastran per questo primo confronto. La rigidezza strutturale e la correzione
aerodinamica rimangono quelle della baseline: si studia il comportamento del
modello corrente, non si certifica l'aereo reale né il contributo del prestress.

Il quesito è separato in due parti:

1. Il protocollo di richiamata genera una componente oscillatoria compatibile
   con il modo BFF, anche senza il rap aggiuntivo?
2. Come cambia la crescita della piccola perturbazione fra 1, 1.3 e 1.6 g?

La seconda domanda usa le nove coppie già eseguite dello sweep `primary` e
la coppia `timestep`. Le nuove run servono a distinguere storia della manovra
e transitorio dovuto al rilascio del SAS. Nessun coefficiente DLM viene ritarato.

## Tre nuove run, un risultato riutilizzato

Primo screening causale a V nominale 67.25 m/s, classe di richiamata 1.6 g.
È il punto più alto dello sweep disponibile, scelto per aumentare
l'osservabilità; non è una ricerca automatica di una conclusione positiva.

| Caso | Comando di manovra | SAS | Origine |
|---|---|---|---|
| `pullup_release` | Picchiata + richiamata originale | Rilascio come baseline | Shadow primary esistente |
| `sham_release` | Ampiezza e pitch-rate comandati nulli | Stesso rilascio | Nuova run |
| `pullup_sas_continuous` | Stessa manovra originale | Sempre attivo | Nuova run |
| `sham_sas_continuous` | Ampiezza e pitch-rate comandati nulli | Sempre attivo | Nuova run |

“Sham” significa controllo senza il comando di manovra, ma con gli stessi
tempi e gli stessi gate degli altri controlli: non è il riferimento a 1 g
precedente con una diversa durata di preparazione. Questo evita di confondere
il tempo di assestamento con l'effetto della manovra.

Tutti i quattro casi hanno rap nullo. Nei casi SAS continuo rimangono attivi
sia i controlli delle superfici sia il damper modale; non si spegne soltanto
uno dei due. Questi casi misurano una risposta a controlli attivi, non un polo
di flutter a controlli disattivati.

I nuovi deck sono copie della shadow realmente eseguita, con sole modifiche
esplicite al comando/gate e all'output. Conservano masse, modi 7–12, smorzamento,
polari, DLM, passo di integrazione 0.01 s e istanti originali della manovra.
L'output parte da t=7 s, prima dell'inizio della picchiata a t=8 s; il tempo
finale è circa 18.15 s. Il passo salvato effettivo si legge dal NetCDF: nella
baseline il meter produce campioni a 0.02 s pur integrando a 0.01 s.

Il confronto comprende la storia picchiata-richiamata: da solo non separa
perfettamente la richiamata dalla precedente picchiata. Se il primo segnale
è convincente, una fase successiva può variare il raccordo del comando e
aggiungere un ingresso senza picchiata. Non è necessario lanciare subito
una matrice più estesa.

## Percorsi

- Input e nuove run: `/home/nicomonzi/ZENO/BFF_PULLUP_CAUSAL_READY`.
- Script: `/home/nicomonzi/X_56/workflows/maneuver_bff`.
- Riferimenti già calcolati: `/mnt/c/Users/Utente/Desktop/BFF_PULLUP_V2`.
- Risultati dell'analisi: `MANOUVER_STIFNESS/results/pullup_bff`.
- Risultati parziali delle due run interrotte su richiesta:
  `/home/nicomonzi/ZENO/BFF_PULLUP_CAUSAL_V1`. Non vengono riutilizzati né cancellati.

I deck contengono percorsi assoluti al modello TESI e il launcher esegue dal
direttorio `BFF_open_loop`. Sono pronti per questo ambiente Linux/WSL; il solo
trasferimento dei `.mbd` su un altro computer non trasferisce le dipendenze.

## Comandi

Per preparare/verificare i tre input senza eseguire alcun solutore:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 run_pullup_causal_controls.py
```

Per lanciare le tre run, al massimo due contemporaneamente:

```bash
cd /home/nicomonzi/X_56/workflows/maneuver_bff
python3 run_pullup_causal_controls.py --execute --jobs 2
```

Attendere tre risultati con `complete: true`. Poi analizzare:

```bash
python3 analyse_pullup_bff.py
```

Il tempo indicativo delle tre run è dell'ordine di 15–30 minuti con due job,
stimato dai precedenti tempi per passo, non misurato su questa nuova campagna.
Dipende da carico CPU e iterazioni. L'analisi richiede molto meno tempo.

`Ctrl+C` chiede l'arresto dei figli e impedisce l'avvio dei casi in coda.
Il launcher rifiuta di sovrascrivere risultati interrotti o input diversi;
per ripartire senza perdere i parziali si usa una nuova cartella esplicita:

```bash
python3 run_pullup_causal_controls.py --output /home/nicomonzi/ZENO/BFF_PULLUP_CAUSAL_RETRY --execute --jobs 2
python3 analyse_pullup_bff.py --controls /home/nicomonzi/ZENO/BFF_PULLUP_CAUSAL_RETRY
```

## Analisi e criteri interpretativi

L'analisi usa direttamente i segnali temporali, con un'oscillazione a frequenza
e crescita stimate più una deriva polinomiale lenta. Non usa soltanto la
pendenza dell'inviluppo Hilbert filtrato su una finestra di pochi cicli.
Controlli indipendenti: uso congiunto di q7 e q7dot, tip simmetriche nel frame
corpo, pitch-rate, ampiezze fra massimi/minimi e sensibilità a finestra/deriva.
La frequenza non è imposta uguale al target DLM.

Per le coppie si analizza `excited-shadow`, che rimuove la manovra comune e
isola la risposta all'impulso. Per le shadow senza rap si separano flessione
lenta ed oscillazione. Un aumento della flessione media non è prova di BFF.

La proiezione delle shadow sulla firma del polo identificato serve a misurare
una componente compatibile con il BFF. Non basta da sola a dimostrare un polo
instabile: occorrono crescita e coerenza del moto rigido/flessibile.
In particolare la proiezione del caso SAS continuo sul polo a SAS spento
non significa che i due sistemi abbiano lo stesso polo.

Le quattro celle permettono di confrontare manovra/non manovra a parità di
gestione del SAS e di verificare l'interazione con il rilascio. Un'oscillazione
presente soltanto dopo il rilascio non si attribuisce automaticamente alla
richiamata. La mancata oscillazione con SAS attivo può invece significare
soppressione da parte del controllo, non assenza della predisposizione al BFF.

L'audit preliminare dei segnali già calcolati ha trovato una forte discordanza
tra il precedente stimatore filtrato e l'adattamento diretto. Perciò i vecchi
onset circa 66.8–66.9 m/s non vanno trattati come verificati. Il nuovo script
conserva il confronto fra metodi, esegue un test sintetico con crescita nota
e non cancella le analisi precedenti. Le conclusioni causali restano in attesa
dei tre controlli; i file di analisi prodotti prima del loro completamento
riportano `causal_controls.status = pending`.

Output attesi:

- `analysis.json`: risultati dettagliati, robustezza, test sintetici e stato controlli;
- `paired_summary.csv`: confronto di crescita/frequenza delle coppie esistenti;
- `existing_trajectories.png`: risposta differenziale e shadow senza rap;
- `causal_controls.png`: confronto dei quattro casi, solo quando completi.

Questa prima campagna non certifica una variazione di flutter fisica dovuta
al prestress e non risolve i limiti della calibrazione DLM. È un test mirato
dell'eccitazione e della crescita del BFF nel modello time-domain disponibile.
