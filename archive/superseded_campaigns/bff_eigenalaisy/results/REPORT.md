# X-56 body-freedom flutter — studio MBDyn 4-DOF rigidi

Densità fissa: `1.039663911 kg/m³` (`9.7284e-08` IPS). Variabile esterna unica: `V_INF`.
Il riferimento NASTRAN separato è `Vf=60.8421 m/s, f=2.0597 Hz` (DLM/PK open-loop); MBDyn è strip-C81 quasi-stazionario, non lineare e closed-loop.
I modi FEM 1–6 sono rigidi e sono esclusi perché il moto rigido è rappresentato dal floating frame MBDyn; il joint esterno rimuove poi soltanto X e Y. I modi flessibili 7–12 coprono 3.217–12.759 Hz.
Il solo joint esterno vincola le traslazioni globali X e Y al CG. Z, roll, pitch e yaw restano liberi; non esistono spinta, throttle o airspeed hold. Rx e Ry sono reazioni ideali e non segnali di controllo.

| V [m/s] | q [Pa] | trim | max abs(Z-Z0) [m] | rigido f [Hz] | sigma rigido [1/s] | BFF f [Hz] | sigma BFF [1/s] | Rx trim [N] | sat. | esito globale |
|---:|---:|:---:|---:|---:|---:|---:|---:|---:|---:|:---|
| 60.8421 | 1924.3 | OK | 0.201 | 0.122 | -0.266 | 1.583 | -0.280 | 124.0 | 0.0% | inconclusive: rigid longitudinal insufficient evidence; BFF stable |

## Controllo e filtri

Le superfici ricevono `delta = delta_trim + delta_SAS + delta_outer + delta_burst + delta_safety`; il burst e' l'unico termine rimosso. Tutti i comandi passano in un attuatore Tustin del primo ordine con tau=0.01 s e limite di correzione ±8 deg.

- scheduling longitudinale: S=1.5 q60/(q_inf+0.5 q60), quindi S=2.0 a 30 m/s, 1.0 a 60 m/s e 0.81 a 70 m/s;
- quota: Kp=4.00e-4 S rad/in, Ki=1.00e-5 S rad/(in s), LP1 a 0.15 Hz; Vz: Kp=1.20e-3 S rad/(in/s), LP1 a 0.35 Hz; entrambi limitati a ±2 deg;
- pitch: Kp=-0.50 S, Ki=-0.008 S; pitch rate: Kp=-0.60 S s; limiti individuali ±6 deg;
- roll: Kp=-0.80, Ki=-0.040; roll rate: Kp=-0.80 s; limiti ±6 deg;
- yaw: Kp=+1.80, Ki=+0.150; yaw rate: Kp=+3.00 s; Vy: Kp=6.0e-4 rad/(in/s), Ki=2.0e-5;
- notch digitale (fs=50 Hz): f0=2.0597 Hz, Q=3, B=[0.95862109,-1.85337926,0.95862109], A=[1,-1.85337926,0.91724218];
- i low-pass longitudinali sono del primo ordine per ridurre il ritardo di fase; il precedente Butterworth a 0.8 Hz resta soltanto sul canale laterale Vy;
- attuatore Tustin: y[k]=0.5 u[k]+0.5 u[k-1];
- safety longitudinale: armato 4 s dopo il burst, deadband |q|=18 deg/s e |pitch|=12 deg, guadagni -0.35 s e -0.25, limite ±4 deg.

## Criteri automatici

Trim valido: X e Y costanti entro 1e-6 in; Rz=0 entro 1e-8 lbf; |roll|<3 deg, |pitch|<8 deg, |alpha|<12 deg, |yaw|<5 deg, |Vz|<2 m/s, errore Vrel<10%, quota<5 m, nessuna saturazione, |Fz| medio <5% del peso e ogni momento medio <75 N m. Il limite pitch ammette l'incidenza fisica necessaria alle basse velocita'; Fx e Fy sono equilibrate esclusivamente da Rx e Ry.
Identificazione valida: run completa, trim valido, almeno 3 cicli, matrix-pencil multi-segnale robusto a nove perturbazioni della finestra, polo nominale interno alla propria CI e saturazione <1%.
Il BFF e' seguito per continuita' in frequenza procedendo dalle velocita' maggiori verso le minori. La classificazione globale e' `stable` soltanto quando sia il modo longitudinale rigido sia il BFF sono identificati e hanno tutta la CI di sigma sotto zero.
