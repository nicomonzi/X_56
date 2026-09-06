import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# =========================================================================
# ANALISI FREQUENZE CON TAGLIO TRANSITORIO E RICERCA MULTI-PICCO
# =========================================================================

NOME_FILE = "/home/nicomonzi/TESI/test/test_mbdyn_aero/output/x.mov" 
DT = 0.001                        
T_CUT = 1  # Secondi da scartare all'inizio della simulazione

# Frequenze estratte dal file f06 di NASTRAN (colonna CYCLES).
# I primi 6 modi sono di corpo rigido (~0 Hz) e vengono ignorati.
FREQUENZE_NASTRAN = [3.217134, 5.302727, 8.705123, 11.16398]

nodi_target = {'990020': 'Tip Sinistra', '991020': 'Tip Destra'}

storia_Z = {'990020': [], '991020': []}
storia_RotY = {'990020': [], '991020': []}

print(f"Lettura file {NOME_FILE}...")

frame_corrente = set()
try:
    with open(NOME_FILE, 'r') as f:
        for riga in f:
            if not riga.strip(): continue
            dati = riga.split()
            nodo_id = dati[0]

            if nodo_id in nodi_target:
                if nodo_id in frame_corrente:
                    frame_corrente.clear()
                frame_corrente.add(nodo_id)

                storia_Z[nodo_id].append(float(dati[3]))
                storia_RotY[nodo_id].append(float(dati[5]))

    n_steps = min(len(storia_Z['990020']), len(storia_Z['991020']))
    tempo = np.arange(n_steps) * DT

    # Calcolo spostamenti netti da t=0 (per il grafico temporale)
    Z_sx = np.array(storia_Z['990020'][:n_steps])
    Z_dx = np.array(storia_Z['991020'][:n_steps])
    delta_Z_sx = Z_sx - Z_sx[0]
    delta_Z_dx = Z_dx - Z_dx[0]

    RotY_sx = np.array(storia_RotY['990020'][:n_steps])
    RotY_dx = np.array(storia_RotY['991020'][:n_steps])
    delta_Pitch_sx = np.degrees(RotY_sx - RotY_sx[0])
    delta_Pitch_dx = np.degrees(RotY_dx - RotY_dx[0])

    # =====================================================================
    # --- TAGLIO DEL TRANSITORIO ---
    # =====================================================================
    idx_cut = np.searchsorted(tempo, T_CUT)
    
    if idx_cut >= n_steps:
        print("Attenzione: La simulazione è più corta del tempo di taglio!")
        idx_cut = 0

    tempo_cut = tempo[idx_cut:]
    
    # Prendiamo solo il segnale a regime e rimuoviamo la sua media statica
    z_sx_regime = delta_Z_sx[idx_cut:] - np.mean(delta_Z_sx[idx_cut:])
    z_dx_regime = delta_Z_dx[idx_cut:] - np.mean(delta_Z_dx[idx_cut:])
    p_sx_regime = delta_Pitch_sx[idx_cut:] - np.mean(delta_Pitch_sx[idx_cut:])
    p_dx_regime = delta_Pitch_dx[idx_cut:] - np.mean(delta_Pitch_dx[idx_cut:])

    # =====================================================================
    # --- ANALISI IN FREQUENZA E RICERCA PICCHI ---
    # =====================================================================
    n_steps_cut = len(tempo_cut)
    freqs = np.fft.rfftfreq(n_steps_cut, d=DT)
    
    # Calcolo distanza minima tra i picchi (es. 0.5 Hz per evitare falsi positivi vicini)
    df = freqs[1] - freqs[0]
    dist_min = max(1, int(0.5 / df))

    def analizza_spettro(segnale):
        fft_vals = np.fft.rfft(segnale)
        mag = np.abs(fft_vals)
        # Troviamo i picchi che siano almeno il 3% del picco massimo
        soglia = np.max(mag) * 0.03
        picchi_idx, _ = find_peaks(mag, height=soglia, distance=dist_min)
        
        # Ordiniamo i picchi in base all'ampiezza (dal più forte al più debole)
        picchi_ordinati = sorted(picchi_idx, key=lambda i: mag[i], reverse=True)
        return mag, picchi_ordinati

    mag_z_sx, picchi_z_sx = analizza_spettro(z_sx_regime)
    mag_z_dx, picchi_z_dx = analizza_spettro(z_dx_regime)
    mag_p_sx, picchi_p_sx = analizza_spettro(p_sx_regime)
    mag_p_dx, picchi_p_dx = analizza_spettro(p_dx_regime)

    # --- FUNZIONE DI SUPPORTO PER LA COMPARAZIONE ---
    def trova_frequenza_nastran_piu_vicina(freq_mbdyn):
        freq_vicina = min(FREQUENZE_NASTRAN, key=lambda x: abs(x - freq_mbdyn))
        errore = abs(freq_mbdyn - freq_vicina) / freq_vicina * 100
        return freq_vicina, errore

    # --- STAMPA DEI RISULTATI ---
    print("\n" + "="*85)
    print(" IDENTIFICAZIONE MODALE E COMPARAZIONE NASTRAN (TRANSITORIO IGNORATO: PRIMI 2.0s)")
    print("="*85)
    
    print(" --- FLESSIONE (Bending) Tip Sinistra ---")
    for i, idx in enumerate(picchi_z_sx[:3]): # Mostra i primi 3 picchi
        f_mbdyn = freqs[idx]
        f_nastran, errore = trova_frequenza_nastran_piu_vicina(f_mbdyn)
        ampiezza_relativa = mag_z_sx[idx]/np.max(mag_z_sx)*100
        print(f" Picco {i+1}: MBDyn = {f_mbdyn:>6.3f} Hz (Amp: {ampiezza_relativa:>5.1f}%) | "
              f"NASTRAN = {f_nastran:>6.3f} Hz | Diff = {errore:>4.1f}%")

    print("\n --- TORSIONE (Pitch) Tip Sinistra ---")
    for i, idx in enumerate(picchi_p_sx[:3]):
        f_mbdyn = freqs[idx]
        f_nastran, errore = trova_frequenza_nastran_piu_vicina(f_mbdyn)
        ampiezza_relativa = mag_p_sx[idx]/np.max(mag_p_sx)*100
        print(f" Picco {i+1}: MBDyn = {f_mbdyn:>6.3f} Hz (Amp: {ampiezza_relativa:>5.1f}%) | "
              f"NASTRAN = {f_nastran:>6.3f} Hz | Diff = {errore:>4.1f}%")
    print("="*85 + "\n")


    # =====================================================================
    # --- GRAFICI ---
    # =====================================================================
    
    # FIGURA 1: Storie Temporali
    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.plot(tempo, delta_Z_sx, 'b-', label='Tip Sinistra', linewidth=1.5)
    ax1.plot(tempo, delta_Z_dx, 'r--', label='Tip Destra', linewidth=1.5, alpha=0.7)
    ax1.axvspan(0, T_CUT, color='gray', alpha=0.3, label='Transitorio Scartato')
    ax1.set_title("Bending (Flessione)", fontweight='bold')
    ax1.set_ylabel(r"$\Delta Z$ [pollici]")
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()

    ax2.plot(tempo, delta_Pitch_sx, 'b-', label='Tip Sinistra', linewidth=1.5)
    ax2.plot(tempo, delta_Pitch_dx, 'r--', label='Tip Destra', linewidth=1.5, alpha=0.7)
    ax2.axvspan(0, T_CUT, color='gray', alpha=0.3)
    ax2.set_title("Torsion (Pitch)", fontweight='bold')
    ax2.set_xlabel("Tempo [s]")
    ax2.set_ylabel(r"$\Delta\theta_y$ [gradi]")
    ax2.grid(True, linestyle='--', alpha=0.7)
    fig1.tight_layout()

    # FIGURA 2: Spettro in Frequenza (FFT)
    fig2, (ax3, ax4) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    ax3.plot(freqs, mag_z_sx, 'b-', label='Spettro Tip Sx')
    ax3.plot(freqs[picchi_z_sx], mag_z_sx[picchi_z_sx], 'rx', markersize=8, label='Modi identificati')
    for idx in picchi_z_sx[:3]:
        ax3.text(freqs[idx]+0.5, mag_z_sx[idx], f"{freqs[idx]:.2f} Hz", fontsize=10)
    # Aggiungo linee verticali per le frequenze di NASTRAN
    for fn in FREQUENZE_NASTRAN:
        if fn < 30:  # Mostra solo quelle nel range del grafico
            ax3.axvline(x=fn, color='green', linestyle=':', alpha=0.6)
    ax3.plot([], [], 'g:', alpha=0.6, label='Modi NASTRAN') # Per la legenda

    ax3.set_title("Spettro Flessione (Segnale a Regime)", fontweight='bold')
    ax3.set_ylabel("Ampiezza Z")
    ax3.set_xlim(0, 30)
    ax3.grid(True, linestyle='--', alpha=0.7)
    ax3.legend()

    ax4.plot(freqs, mag_p_sx, 'b-', label='Spettro Tip Sx')
    ax4.plot(freqs[picchi_p_sx], mag_p_sx[picchi_p_sx], 'rx', markersize=8)
    for idx in picchi_p_sx[:3]:
        ax4.text(freqs[idx]+0.5, mag_p_sx[idx], f"{freqs[idx]:.2f} Hz", fontsize=10)
    for fn in FREQUENZE_NASTRAN:
        if fn < 30:
            ax4.axvline(x=fn, color='green', linestyle=':', alpha=0.6)

    ax4.set_title("Spettro Torsione (Segnale a Regime)", fontweight='bold')
    ax4.set_xlabel("Frequenza [Hz]")
    ax4.set_ylabel(r"Ampiezza $\theta_y$")
    ax4.set_xlim(0, 30) 
    ax4.grid(True, linestyle='--', alpha=0.7)

    fig2.tight_layout()
    plt.show()

except FileNotFoundError:
    print(f"Errore: Impossibile trovare il file '{NOME_FILE}'.")