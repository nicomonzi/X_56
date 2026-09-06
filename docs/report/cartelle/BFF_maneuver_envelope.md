# Rapporto: `BFF_maneuver_envelope`

## Ruolo e stato

Contiene il generatore parametrico e gli analizzatori della campagna BFF durante
dive–pull-up. È una **dipendenza attiva** di `MANOUVER_STIFNESS`: non va
archiviata come semplice prototipo. Lo snapshot ha 21 file e circa 148 MB,
quasi tutti concentrati in due NetCDF da 77,2 MB dentro
`output_n1_verification`.

La working tree è non pulita: quattro file tracciati sono modificati, tre vecchi
script risultano rimossi e i nuovi workflow dive–pull-up/output sono non
tracciati. È quindi una zona ad alto rischio di perdita o commistione.

## Codice

- `maneuver_case.py` definisce `ManeuverPoint`, carica la baseline
  `BFF_open_loop`, risolve i tempi della manovra, sostituisce setpoint e costanti,
  renderizza l'input MBDyn ed esegue un singolo caso;
- `run_dive_pullup_sweep.py` costruisce le coppie shadow/excited, scrive il
  manifest e non esegue nulla senza `--execute`;
- `analyse_dive_pullup.py` misura traiettoria, onset, coerenza tip/modo e crea
  grafici riepilogativi;
- `analyse_time_domain_pairs.py` ricostruisce tip simmetrica, leakage delle
  superfici e metriche di volo di coppie steady/manovra;
- `compare_paired_response.py` sottrae i segnali modali delle coppie e stima
  crescita/frequenza;
- `study_config.json` contiene griglia, tempi e criteri.

Il disegno usa una shadow e una excited per lo stesso comando cinematico. La
sola differenza è il piccolo RAP simmetrico WF4; `excited-shadow` elimina al
primo ordine la risposta forzata comune. Durante SAS-off tutte le superfici sono
in sample-and-hold. Il carico `n` non è un setpoint: viene ricostruito a
posteriori dalla traiettoria.

## Avanzamento

Il README descrive 24 traiettorie su 65/66/67 m/s e quattro aggressività, ma la
linea più recente ha ristretto e corretto il piano dentro
`MANOUVER_STIFNESS`. I due NetCDF locali a 66,75 m/s e n=1 sono una verifica,
non l'intera campagna. Gli output principali citati dal codice vivono sul
Desktop Windows e non sono inventariati qui.

## Obsolescenza e azioni

La cartella non è obsoleta, ma le funzioni esecutive e scientifiche vanno
separate: mantenere qui il generatore di manovra riusabile; spostare matrici,
criteri e risultati nella campagna che li usa. Consolidare subito lo stato Git,
decidere esplicitamente il destino dei tre script rimossi e non versionare i
NetCDF grezzi senza una politica per grandi file.

