# Dizionario dei dati

I CSV sono UTF-8 con intestazione e vengono rigenerati da
`../aggiorna_report.py`. Le dimensioni sono in byte; le date sono UTC ISO 8601.

- `inventario_file.csv`: una riga per ogni file regolare, esclusi `.git` e
  `report`; `tracciato_git=no` comprende sia file ignorati sia non tracciati.
- `duplicati_sha256.csv`: una riga per ogni occorrenza duplicata. Le righe con
  lo stesso `gruppo` hanno hash e contenuto identici. Il campo
  `byte_ridondanti_massimi` è ripetuto sulle righe del gruppo e non deve essere
  sommato senza prima raggruppare.
- `riepilogo_cartelle.csv`: aggregazione sulle 17 cartelle di primo livello.
- `inventario_sottocartelle.csv`: una riga per ogni directory annidata,
  incluse directory vuote; i conteggi ricorsivi includono tutti i discendenti.
- `catalogo_codice_python.csv`: inventario AST degli script di progetto;
  ambienti virtuali, build, site-packages e cache sono esclusi.
- `cronologia_git.csv`: commit ordinati dal più recente al più vecchio.
- `snapshot.json`: hash e metadati usati solo per il confronto incrementale.

Il massimo recuperabile dai duplicati non considera la necessità di copie
locali per include, casi autonomi o riproducibilità. Va usato come indicatore
per una revisione, non come elenco di cancellazione.

