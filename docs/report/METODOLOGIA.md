# Metodologia dell'audit

## Perimetro

La scansione copre tutti i file regolari sotto `/home/nicomonzi/TESI`, escluse
la directory interna `.git` e la directory `report`. I link simbolici non sono
contati come file regolari; i tre link individuati sono riportati nel riepilogo
di snapshot. Le directory `.agents` e `.codex` sono vuote e quindi non hanno un
rapporto separato. Il contenuto di `/home/nicomonzi/ZENO` e dei percorsi Windows
citati dai launcher è stato usato solo quando già descritto nei documenti del
progetto; non fa parte dell'inventario quantitativo di TESI.

## Fonti

Le conclusioni combinano cinque fonti, in quest'ordine di affidabilità:

1. file correnti e risultati numerici presenti su disco;
2. README, audit, manifest e report generati dagli script del progetto;
3. codice Python, shell, XML preCICE e input MBDyn/Nastran/DUST;
4. stato e cronologia Git;
5. cronologia operativa disponibile nella conversazione, usata come contesto e
   non come sostituto dell'evidenza su disco.

Quando un README storico è contraddetto da una revisione più recente, prevale
la revisione. È il caso di `MANOUVER_STIFNESS`: la rettifica del 5 settembre
invalida il vecchio gate sul prestress e il nuovo stimatore mette in attesa la
conclusione sull'onset da manovra.

## Definizione di duplicato

Due file sono duplicati solo se hanno dimensione maggiore di zero e lo stesso
SHA-256, cioè sono identici byte per byte. Sono fornite tre quantità diverse:

- **gruppi duplicati**: contenuti distinti che compaiono almeno due volte;
- **occorrenze**: tutti i file appartenenti a quei gruppi;
- **copie extra**: ogni occorrenza oltre la prima; è il numero più vicino a
  “quanti file sono ripetuti”.

I 252 file vuoti non sono inclusi: considerarli un unico contenuto duplicato
aggiungerebbe 251 copie ma non rappresenterebbe spazio né informazione
replicata. Il “massimo recuperabile” presume di conservare una sola copia per
hash; non è una raccomandazione di cancellazione. Molte copie FEM, C81, BULK e
input locali sono intenzionali per la riproducibilità autonoma dei casi.

## Criterio di obsolescenza

`Obsoleta` non significa “senza valore scientifico”. Le etichette usate sono:

- **attiva/canonica**: è una sorgente corrente o una dipendenza dei workflow
  recenti;
- **attiva con blocchi**: il codice è corrente, ma una conclusione o produzione
  è esplicitamente sospesa;
- **supporto/validazione**: non è il workflow finale, ma produce una baseline,
  una verifica o un asset ancora usato;
- **archivio/superseded**: è stato sostituito per l'esecuzione corrente, ma va
  conservato finché la provenienza è migrata e verificata;
- **generata/rigenerabile**: cache, build o output ricostruibile; è il candidato
  più sicuro alla pulizia, dopo verifica del manifest e backup.

La classificazione è basata su dipendenze nel codice, README, data e natura dei
risultati. Nessuna directory è stata spostata o cancellata durante l'audit.

## Limiti

- Non sono state rilanciate simulazioni MBDyn, DUST o Nastran: i risultati sono
  valutati tramite file, log e report esistenti.
- Il repository è deliberatamente fotografato in stato non pulito. La nuova
  cartella `MANOUVER_STIFNESS` è ancora fuori da Git e il ramo
  `BFF_maneuver_envelope` contiene modifiche, rimozioni e nuovi output.
- I messaggi Git sono spesso generici (`ok`, `update`, `m`); la cronologia
  tecnica è quindi ricostruita soprattutto dal contenuto e dai documenti.
- Le date `mtime` possono cambiare per copia o checkout; non sono usate da sole
  per decidere l'avanzamento.

