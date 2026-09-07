import matplotlib.pyplot as plt
import numpy as np

# =========================================================================
# ANALISI STORIA TEMPORALE DEFORMAZIONI X-56A (MBDyn)
# =========================================================================

NOME_FILE = "/home/nicomonzi/TESI/test/test_mbdyn_aero/output/x.mov" # Assicurati che il nome sia corretto
DT = 0.001                         # Il time-step che hai impostato in MBDyn

# Nodi da tracciare (le due Tip prima delle winglet)
nodi_target = {'990020': 'Tip Sinistra', '991020': 'Tip Destra'}

# Dizionari per salvare la storia temporale
storia_Z = {'990020': [], '991020': []}
storia_RotY = {'990020': [], '991020': []}

print(f"Estrazione dati per i nodi {list(nodi_target.keys())} dal file {NOME_FILE}...")

frame_corrente = set()
step_count = 0

try:
    with open(NOME_FILE, 'r') as f:
        for riga in f:
            if not riga.strip(): continue
            dati = riga.split()
            nodo_id = dati[0]

            if nodo_id in nodi_target:
                # Logica per separare i time-step: 
                # se abbiamo già visto questo nodo, è iniziato un nuovo fotogramma
                if nodo_id in frame_corrente:
                    frame_corrente.clear()
                    step_count += 1
                
                frame_corrente.add(nodo_id)

                # Estrazione dati (gli indici Python partono da 0)
                # Dati[3] = Z (Flessione)
                # Dati[5] = RotY (Componente Y del vettore orientazione -> Torsione)
                Z = float(dati[3])
                RotY = float(dati[5])

                storia_Z[nodo_id].append(Z)
                storia_RotY[nodo_id].append(RotY)

    # Assicuriamoci che i vettori abbiano la stessa lunghezza
    n_steps = min(len(storia_Z['990020']), len(storia_Z['991020']))
    tempo = np.arange(n_steps) * DT

    print(f"Estratti {n_steps} passi temporali (Fino a t = {tempo[-1]:.3f} s).")

    # --- CALCOLO DEFORMAZIONI RELATIVE (Delta da t=0) ---
    Z_sx = np.array(storia_Z['990020'][:n_steps])
    Z_dx = np.array(storia_Z['991020'][:n_steps])
    
    # Sottraiamo il valore a t=0 per avere la flessione netta in pollici
    delta_Z_sx = Z_sx - Z_sx[0]
    delta_Z_dx = Z_dx - Z_dx[0]

    RotY_sx = np.array(storia_RotY['990020'][:n_steps])
    RotY_dx = np.array(storia_RotY['991020'][:n_steps])
    
    # Sottraiamo il valore a t=0 per avere la torsione netta e convertiamo in GRADI
    delta_Pitch_sx = np.degrees(RotY_sx - RotY_sx[0])
    delta_Pitch_dx = np.degrees(RotY_dx - RotY_dx[0])

    # --- CREAZIONE DEI GRAFICI ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # GRAFICO 1: FLESSIONE (BENDING)
    ax1.plot(tempo, delta_Z_sx, 'b-', label='Tip Sinistra (990020)', linewidth=2)
    ax1.plot(tempo, delta_Z_dx, 'r--', label='Tip Destra (991020)', linewidth=2)
    ax1.set_title("Bending", fontsize=14, fontweight='bold')
    ax1.set_ylabel("$\Delta Z$ [pollici]", fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend(loc='upper right')

    # GRAFICO 2: TORSIONE (PITCH)
    ax2.plot(tempo, delta_Pitch_sx, 'b-', label='Tip Sinistra (990020)', linewidth=2)
    ax2.plot(tempo, delta_Pitch_dx, 'r--', label='Tip Destra (991020)', linewidth=2)
    ax2.set_title("Torsion", fontsize=14, fontweight='bold')
    ax2.set_xlabel("Tempo [s]", fontsize=12)
    ax2.set_ylabel("$\Delta\\theta_y$ [gradi]", fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

except FileNotFoundError:
    print(f"Errore: Impossibile trovare il file '{NOME_FILE}'.")
except IndexError:
    print("Errore formato: Il file .mov non sembra avere le colonne dell'orientazione.")