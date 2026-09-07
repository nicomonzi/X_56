import os
import re
import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.signal import find_peaks

# =============================================================================
# SETUP DEI PARAMETRI
# =============================================================================
template_file = 'main_x56.mbd'   # Il tuo file MBDyn originale
output_dir = 'output'                # Cartella dove finiranno i risultati
velocities_ms = np.arange(55, 60, 1) # Range di velocità da 34 a 44 m/s

# Nodi di interesse
node_center = '990001'
node_tip_R = '990020'
node_tip_L = '991020'

# Inizializzazione vettori per i diagrammi riassuntivi
dampings = []
frequencies = []
velocities_plot = []

# Crea la cartella di output se non esiste
os.makedirs(output_dir, exist_ok=True)

# =============================================================================
# FUNZIONI DI SUPPORTO
# =============================================================================
def calc_damping_frequency(time_array, signal):
    """
    Calcola smorzamento (g) e frequenza (f) usando il decremento logaritmico 
    sui picchi del segnale nel dominio del tempo.
    """
    # Trova i picchi nel segnale
    peaks, _ = find_peaks(signal)
    
    if len(peaks) < 3:
        return np.nan, np.nan # Non ci sono abbastanza cicli per una stima affidabile
        
    t_peaks = time_array[peaks]
    y_peaks = signal[peaks]
    
    # Frequenza: media dell'inverso del periodo tra i picchi
    periods = np.diff(t_peaks)
    f = 1.0 / np.mean(periods)
    
    # Decremento logaritmico delta = (1 / n) * ln(x_0 / x_n)
    n = len(y_peaks) - 1
    
    # Shift del segnale se non è centrato in zero (per evitare log negativi)
    if y_peaks[-1] <= 0 or y_peaks[0] <= 0:
        y_peaks = y_peaks - np.mean(signal)
        
    # Calcolo delta prendendo il valore assoluto per sicurezza
    delta = (1 / n) * np.log(np.abs(y_peaks[0] / y_peaks[-1]))
    
    # Smorzamento adimensionale zeta e conversione in g (g = 2*zeta)
    zeta = delta / np.sqrt(4 * np.pi**2 + delta**2)
    g = 2 * zeta 
    
    return g, f

# =============================================================================
# LOOP DI SIMULAZIONE E ANALISI
# =============================================================================
# Leggi il template una volta sola in memoria
if not os.path.exists(template_file):
    raise FileNotFoundError(f"Errore: Il file {template_file} non e' stato trovato nella directory corrente.")

with open(template_file, 'r') as f:
    template_content = f.read()

for V in velocities_ms:
    print(f"\n--- Avvio simulazione a V = {V} m/s ---")
    start_time = time.time()
    
    # 1. Aggiorna la velocità nel file MBDyn (mantiene la conversione * m2in)
    new_content = re.sub(
        r'set: const real VINF\s*=\s*[\d\.\*m2in\s]+;', 
        f'set: const real VINF = {V} * m2in;', 
        template_content
    )
    
    run_file_v = f'x56_run_V{V}.mbd'
    with open(run_file_v, 'w') as f:
        f.write(new_content)
        
    # 2. Lancia MBDyn
    out_prefix = f'out_V{V}'
    out_path = os.path.join(output_dir, out_prefix)
    comando_mbdyn = ['mbdyn', '-f', run_file_v, '-o', out_path]
    
    print(f"Esecuzione in corso: {' '.join(comando_mbdyn)}")
    subprocess.run(comando_mbdyn, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    elapsed_time = time.time() - start_time
    print(f"Simulazione completata in {elapsed_time:.1f} secondi.")
    
    # 3. Leggi i risultati da NetCDF
    nc_file = f'{out_path}.nc'
    if not os.path.exists(nc_file):
        print(f"Errore: file {nc_file} non trovato. MBDyn potrebbe aver fallito la run.")
        continue
        
    try:
        ds = Dataset(nc_file, 'r')
        time_array = ds.variables['time'][:]
        
        # Estrazione Z (indice 2) e Pitch (indice 1)
        tip_z = ds.variables[f'node.struct.{node_tip_R}.X'][:, 2] 
        tip_pitch = ds.variables[f'node.struct.{node_tip_R}.Phi'][:, 1] 
        
        center_z = ds.variables[f'node.struct.{node_center}.X'][:, 2]
        center_pitch = ds.variables[f'node.struct.{node_center}.Phi'][:, 1]
        
    except KeyError as e:
        print(f"Errore nella lettura del NetCDF per V={V}. Verifica i nomi delle variabili: {e}")
        ds.close()
        continue 
        
    # 4. Calcola smorzamento e frequenza (tagliando il primo 10% di transitorio numerico)
    transient_cut = int(len(time_array) * 0.1) 
    g, f_val = calc_damping_frequency(time_array[transient_cut:], tip_z[transient_cut:])
    
    dampings.append(g)
    frequencies.append(f_val)
    velocities_plot.append(V)
    
    # 5. Plot Time-History per la singola velocità
    plt.figure(figsize=(12, 5))
    
    # Plot Deflessione Z
    plt.subplot(1, 2, 1)
    plt.plot(time_array, tip_z, label=f'Tip R Z-disp', color='blue')
    plt.plot(time_array, center_z, label='Center Z-disp', color='gray', alpha=0.5)
    plt.xlabel('Time [s]')
    plt.ylabel('Deflection Z [in]')
    plt.title(f'Wing Deflection (V={V} m/s)')
    plt.legend()
    plt.grid(True)
    
    # Plot Rotazione Y (Pitch)
    plt.subplot(1, 2, 2)
    plt.plot(time_array, tip_pitch, label=f'Tip R Pitch', color='red')
    plt.plot(time_array, center_pitch, label='Center Pitch', color='gray', alpha=0.5)
    plt.xlabel('Time [s]')
    plt.ylabel('Rotation Y [rad]')
    plt.title(f'Wing Rotation (V={V} m/s)')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    # Salva il grafico nella cartella output
    plt.savefig(os.path.join(output_dir, f'TimeHistory_V{V}.png'))
    plt.close() # Importante per liberare memoria durante il loop
    
    # Chiudi il dataset NetCDF
    ds.close()
    
    # Pulizia: elimina il file di input temporaneo
    os.remove(run_file_v)

# =============================================================================
# PLOT DIAGRAMMI RIASSUNTIVI V-g e V-f
# =============================================================================
if len(velocities_plot) > 0:
    plt.figure(figsize=(14, 6))

    # V-g Plot (Smorzamento)
    plt.subplot(1, 2, 1)
    plt.plot(velocities_plot, dampings, 'o-', color='crimson', linewidth=2)
    plt.axhline(0, color='black', linestyle='--') # Linea di instabilità (Flutter)
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Damping (g)')
    plt.title('V-g Diagram (Tip Dominant Mode)')
    plt.grid(True)

    # V-f Plot (Frequenza)
    plt.subplot(1, 2, 2)
    plt.plot(velocities_plot, frequencies, 's-', color='navy', linewidth=2)
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Frequency [Hz]')
    plt.title('V-f Diagram (Tip Dominant Mode)')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Vg_Vf_Diagrams.png'))
    plt.show()

    print("\nAnalisi completata con successo!")
    print(f"Tutti i grafici e i file di output MBDyn sono stati salvati nella cartella '{output_dir}'.")
else:
    print("\nNessun dato valido estratto per plottare i diagrammi V-g e V-f.")