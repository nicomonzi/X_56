# Analisi trim accoppiato a 55 m/s

La simulazione è completa fino a 7.996 s ed è numericamente stabile,
ma **il trim non è stato raggiunto**.

- Pitch congelato: 0.8790 deg.
- BFL/BFR congelati: 10.0000 / 10.0000 deg;
  saturazione a +10 deg da t=5.672 s.
- Fz aerodinamica media 6--8 s: 683.262 lbf,
  contro 419.440 lbf.
- Reazione verticale media residua: 263.822 lbf.
- Momento aerodinamico medio al CG: 1910.191 lbf in.
- Reazione di beccheggio media residua: 1910.645 lbf in.
- Deflessione simmetrica media delle estremità: 0.6018 in.

Le pendenze nella finestra congelata sono piccole: il sistema si è assestato, ma
su un equilibrio imposto dai vincoli e non su una condizione di volo trimmata.
Il comando dei body flap procede nel verso che aumenta la deflessione positiva
fino alla saturazione senza annullare My: prima di una nuova run va corretto il
segno/il Jacobiano del canale BFL/BFR e va resa coerente la condizione iniziale.
Il seed usato (0.770 deg, 3.049 deg circa) non è coerente
con il trim DUST precedentemente trovato vicino a 0.527 deg e -3.301 deg.

MBDyn converge con una media di 1.999
iterazioni e un massimo di 2. Dopo il
congelamento non si osserva crescita delle coordinate modali; il modo 7 domina
la deformazione statica. Questo indica regolarità numerica del trim vincolato,
non costituisce una verifica di flutter o del SAS.
