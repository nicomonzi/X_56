# Filtri, notch, attuatori e taratura del controllo

## Convenzione MBDyn

Il `discrete filter drive` di MBDyn implementa

\[
y_k=\sum_{i=1}^{n_a}a_i y_{k-i}+\sum_{i=0}^{n_b}b_i u_{k-i}.
\]

La documentazione locale è in
`/home/nicomonzi/src/mbdyn/manual/input/general.tex`, righe 2446–2497. La
nota fondamentale è che il filtro richiede un passo temporale fisso. Qui
`dt=0.02 s`, quindi la frequenza di campionamento è `fs=50 Hz` e Nyquist è
`25 Hz`.

Il test ufficiale MBDyn è
`/home/nicomonzi/src/mbdyn/tests/drives/drivecallers/discretefilter`, righe
66–73: applica un Butterworth del secondo ordine a una forza sinusoidale.

## Notch

Un notch è un filtro elimina-banda molto stretto: attenua fortemente una
frequenza scelta lasciando quasi invariato il resto. Non stabilizza da solo il
velivolo; riduce quanto il SAS reagisce a una particolare oscillazione.

Nel modello il notch è centrato su `f0=2.0597 Hz`, la frequenza BFF NASTRAN,
con `Q=3`. La larghezza indicativa è `f0/Q = 0.687 Hz`. I coefficienti sono
stati ottenuti con SciPy:

```python
from scipy.signal import iirnotch
B, A = iirnotch(2.0597, Q=3.0, fs=50.0)
a1, a2 = -A[1], -A[2]
b0, b1, b2 = B
```

Risultato MBDyn:

```text
a1 =  1.853379259531933
a2 = -0.917242175082742
b0 =  0.958621087541371
b1 = -1.853379259531933
b2 =  0.958621087541371
```

Il notch è applicato ai ratei `p`, `q`, `r`, non agli angoli. Il ramo MBDyn
osservato attorno a 1.48–1.60 Hz non coincide col centro del notch: perciò il
SAS continua a partecipare a quel ramo. Questa scelta segue la richiesta di
rendere il SAS trasparente soprattutto alla frequenza BFF NASTRAN.

## Low-pass dei loop lenti

Il precedente secondo ordine a 0.8 Hz cadeva nella stessa banda del modo
longitudinale rigido osservato e introduceva troppo ritardo di fase. Quota e
`Vz` usano ora due Butterworth del primo ordine distinti:

- quota: `fc=0.15 Hz`, per correggere soltanto la deriva lenta;
- `Vz`: `fc=0.35 Hz`, per conservare damping sul moto rigido e attenuare di
  oltre 12 dB il segnale nella banda BFF 1.5–1.7 Hz.

I coefficienti derivano da:

```python
from scipy.signal import butter
B, A = butter(1, fc, btype="low", fs=50.0)
a1 = -A[1]
b0, b1 = B
```

Risultati MBDyn:

```text
quota 0.15 Hz: a1=0.981325890492688, b0=b1=0.009337054753656
Vz    0.35 Hz: a1=0.956957321922681, b0=b1=0.021521339038659
```

Il Butterworth del secondo ordine a 0.8 Hz rimane soltanto sul canale
laterale lento `Vy`.

## Attuatore del primo ordine

L’attuatore rappresenta il ritardo tra comando e deflessione, invece di far
muovere istantaneamente la superficie. Il modello continuo è

\[
\tau \dot\delta+\delta=\delta_c,\qquad \tau=0.01\;s.
\]

Con Tustin e `dt=0.02 s`:

\[
a_1=\frac{2\tau-dt}{2\tau+dt}=0,
\quad b_0=b_1=\frac{dt}{2\tau+dt}=0.5.
\]

In MBDyn diventa quindi

```text
discrete filter, 1, 0.0, 0.5, 1, 0.5, <comando>
```

ossia `delta[k]=0.5*command[k]+0.5*command[k-1]`. È un primo ordine
continuo discretizzato, non un attuatore idraulico MBDyn. L’attuatore
idraulico presente in `manual/input/elemhydr.tex` è un elemento diverso e non
serve per questa simulazione aeroelastica.

## Come sono stati scelti i guadagni quota

La taratura è stata fatta conservando separazione di banda e autorità fisica
delle sole superfici:

1. il seed di trim pitch/body-flap è stato corretto in funzione di
   `TRIM_SCALE`, eliminando il comando statico di circa 3.6° che il SAS doveva
   fornire a 30 m/s;
2. l'outer loop è stato portato sotto la banda dei modi dinamici con i due
   filtri del primo ordine descritti sopra;
3. i guadagni longitudinali sono schedulati mediante
   `S=1.5*q60/(q_inf+0.5*q60)`: `S=2` a 30 m/s, `S=1` a 60 m/s e circa
   `S=0.81` a 70 m/s;
4. la taratura è `Kp_h=4.00e-4*S rad/in`,
   `Ki_h=1.00e-5*S rad/(in s)`, `Kp_Vz=1.20e-3*S rad/(in/s)`;
5. anche pitch e pitch-rate sono schedulati:
   `Kp_theta=-0.50*S`, `Ki_theta=-0.008*S`, `Kp_q=-0.60*S s`.

Nelle verifiche 30–62.5 m/s la quota massima è rimasta entro circa 0.30 m,
con errore finale 0.05–0.12 m, superfici sotto circa 1.2° e nessuna
saturazione. A 65 m/s il BFF porta invece alla saturazione: il controllo quota
non viene usato come active-flutter-suppression e non deve mascherare questa
frontiera fisica.

## PID e documentazione locale

- Sintassi e anti-windup del modulo PID:
  `/home/nicomonzi/src/mbdyn/manual/input/module-pid.tex`, righe 34–71.
- Esempio moderno di altitude PID:
  `/home/nicomonzi/src/mbdyn/modules/module-pid/rigid_rotor_pid.mbd`, righe
  232–246.
- Esempio PID storico tramite state-space:
  `/home/nicomonzi/src/mbdyn/tests/genel/statespace/pid`.
- Sorgente del modulo:
  `/home/nicomonzi/src/mbdyn/modules/module-pid/module-pid.cc`.
- Manuale input in formato sorgente:
  `/home/nicomonzi/src/mbdyn/manual/input/`.
- Per generare il PDF, il target è descritto in
  `/home/nicomonzi/src/mbdyn/manual/Makefile.am`, righe 73–90, e richiama
  `./make.sh --tgt input` dalla directory `manual`.
