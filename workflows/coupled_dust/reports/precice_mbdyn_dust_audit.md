# Audit accoppiamento MBDyn - preCICE - DUST

Data: 2026-09-08

## Riferimento

Il controllo segue il capitolo 4, pagine 41-46, del manuale DUST fornito
(`DUST_user_manual.pdf`): componente `coupled`, accoppiamento `rbf`, cinematica
da MBDyn a DUST, carichi da DUST a MBDyn, mapping `consistent` in andata e
`conservative` al ritorno, schema seriale implicito e accelerazione Aitken.

Il manuale e del 2023 e mostra la sintassi XML preCICE precedente. La preCICE
3.2.0 installata non accetta l'attributo `valid-digits` in
`time-window-size`; l'attributo e stato quindi rimosso dal solo file XML v3.
Il file legacy v2 non e stato modificato.

## Toolchain usata

- MBDyn: `/usr/local/mbdyn/bin/mbdyn`
- interfaccia socket MBDyn: `/usr/local/mbdyn/libexec/mbpy/mbc_py_interface`
- DUST patched: `/home/monzani/dust-patched/build-user/bin`
- libreria preCICE caricata da DUST e dal binding Python: `libprecice.so.3.2.0`
- API Python: `Participant` v3

I percorsi sono registrati nella configurazione locale ignorata da Git
`config/machine.env`.

## Topologia e sistemi di riferimento

| controllo | risultato |
|---|---:|
| nodi nella external force MBDyn | 59, tutti unici |
| nodi strutturali principali | 39 |
| nodi estremi hinge | 20, due per ciascuna delle 10 superfici |
| `refConfigNodes.in` contro `coupling_nodes.in` | differenza massima 0 |
| nodi nel file HDF5 FINE contro i due file ASCII | differenza massima 0 |
| mesh FINE | 7320 pannelli, 7503 punti |

La geometria aerodinamica, i nodi MBDyn e i nodi di coupling sono tutti gia
espressi nel sistema globale IPS. Per questo caso la matrice di orientazione
identita e corretta; non serve la trasformazione tra assi locali mostrata
nell'esempio ad ala separata del manuale.

## Scambio dati

- MBDyn scrive `Position`, `Rotation`, `Velocity`, `AngularVelocity`.
- DUST scrive `Force`, `Moment`.
- Il trasferimento cinematico e `nearest-neighbor consistent`.
- Il trasferimento dei carichi e `nearest-neighbor conservative`.
- Il socket MBDyn usa `coupling, tight`, `sorted, yes`, orientation vector e
  accelerazioni, coerentemente con l'adapter.

La directory di handshake preCICE e stata resa assoluta e specifica per il
caso. I soli metadati di connessione residui da una run interrotta vengono
rimossi prima dell'avvio.

## Prove eseguite

### Parsing

MBDyn e DUST raggiungono entrambi l'attesa di coupling con:

- COARSE: 2440 pannelli, 2583 punti, 59 nodi;
- FINE: 7320 pannelli, 7503 punti, 59 nodi e 10 hinge.

### Smoke COARSE e FINE

Entrambi i test avanzano da 0 a 0.1 s con passo 0.002 s. Il test FINE ha
prodotto 50 stati strutturali accettati in 906.83 s:

- 3 o 4 iterazioni implicite per finestra, massimo configurato 20;
- tutti i campi cinematici e i carichi finiti;
- DUST termina con `Computations Finished`;
- MBDyn e DUST terminano con codice zero;
- output HDF5, CSV e VTU generati correttamente.

Gli avvisi Aitken alla prima finestra sono dovuti a sottovettori cinematici
inizialmente invariati. Non si traducono in divergenza o valori non finiti.
Il messaggio `got ABORT from peer` compare durante la chiusura MPI dopo la
terminazione regolare dei due solver e non modifica i codici di uscita.

### Conservazione dei carichi

Sono stati sommati direttamente i campi `dF` e `dMom` del risultato HDF5
DUST a 0.080 s e confrontati con i carichi nodali ricevuti da MBDyn nella
finestra successiva a 0.082 s, come richiesto dall'ordine seriale dello schema.

| grandezza | errore relativo |
|---|---:|
| forza risultante | 4.35e-8 % |
| momento risultante all'origine | 6.71e-4 % |

La seconda differenza include l'approssimazione usata nell'audit per
ricostruire il centro geometrico dei pannelli. Il trasferimento conserva
quindi risultante, momento e segni.

## Esito

L'adapter preCICE e correttamente collegato al modello MBDyn e alla mesh DUST
FINE da 30 elementi in corda. L'accoppiamento strutturale e idoneo per passare
alla preparazione dei casi open-loop.

Resta da eseguire come prova distinta il comando non nullo delle dieci
superfici mobili. Gli smoke test qui documentati mantengono le hinge neutre e
validano la loro topologia, ma non ancora la cinematica a deflessione finita.
