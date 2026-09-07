import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import sys

# ==============================================================================
# CONFIGURAZIONE
# ==============================================================================
FILE_NC = "/home/nicomonzi/TESI/PID_test/output/x.nc"  
ID_BASE = "990000"
ID_ALA = "990001"
ID_PID = "20"
SETPOINT_ROLLIO = 10.0  # Gradi

# ==============================================================================
# LETTURA DEL FILE NetCDF
# ==============================================================================
try:
    dataset = nc.Dataset(FILE_NC, 'r')
    print(f"File '{FILE_NC}' caricato con successo.\n")
except FileNotFoundError:
    print(f"ERRORE: Impossibile trovare '{FILE_NC}'.")
    sys.exit()

try:
    time = dataset.variables['time'][:]
    num_steps = len(time)
except KeyError:
    print("ERRORE: Variabile 'time' non trovata.")
    sys.exit()

# --- FUNZIONI DI ESTRAZIONE SICURA ---
def estrai_dati(var_nc, steps):
    dati = var_nc[:] 
    if len(dati.shape) > 1:
        return dati[:, 0]
    else:
        if len(dati) == steps * 3:
            return dati[0::3]
        elif len(dati) == steps:
            return dati
        else:
            return dati[:steps]

def estrai_angoli_3d(dataset, id_nodo, steps):
    var_name_E = f'node.struct.{id_nodo}.E'
    var_name_Phi = f'node.struct.{id_nodo}.Phi'
    
    var_name = var_name_E if var_name_E in dataset.variables else None
    if not var_name:
        var_name = var_name_Phi if var_name_Phi in dataset.variables else None
        
    if var_name:
        dati = dataset.variables[var_name][:]
        if len(dati.shape) == 2 and dati.shape[1] >= 3:
            e1, e2, e3 = dati[:, 0], dati[:, 1], dati[:, 2]
        elif len(dati.shape) == 1 and len(dati) == steps * 3:
            e1, e2, e3 = dati[0::3], dati[1::3], dati[2::3]
        else:
            return np.zeros(steps), np.zeros(steps), np.zeros(steps)
        return (np.degrees(np.asarray(e1, dtype=float)),
                np.degrees(np.asarray(e2, dtype=float)),
                np.degrees(np.asarray(e3, dtype=float)))
    else:
        return np.zeros(steps), np.zeros(steps), np.zeros(steps)

# 1. Estrazione Angoli di Eulero (Base e Ala)
base_E1, base_E2, base_E3 = estrai_angoli_3d(dataset, ID_BASE, num_steps)
ala_E1, ala_E2, ala_E3 = estrai_angoli_3d(dataset, ID_ALA, num_steps)

# 2. Estrazione Velocità di Rollio (Base)
base_omega1 = np.zeros(num_steps)
if f'node.struct.{ID_BASE}.Omega' in dataset.variables:
    base_omega1_rad = estrai_dati(dataset.variables[f'node.struct.{ID_BASE}.Omega'], num_steps)
    base_omega1 = np.degrees(np.asarray(base_omega1_rad, dtype=float))

# 3. Estrazione Azione del PID (Flap)
flap_cmd_rad = np.zeros(num_steps)
nome_var_pid = None
for var in dataset.variables.keys():
    if f'.{ID_PID}.' in var and np.issubdtype(dataset.variables[var].dtype, np.number):
        nome_var_pid = var
        break

if nome_var_pid:
    flap_cmd_rad = estrai_dati(dataset.variables[nome_var_pid], num_steps)
    
flap_cmd_deg = np.degrees(np.asarray(flap_cmd_rad, dtype=float))

dataset.close()

# ==============================================================================
# STAMPA DI DIAGNOSTICA SUL TERMINALE
# ==============================================================================
print("\n" + "="*70)
print(" 📊 RIASSUNTO DATI ESTRATTI (TELEMETRIA RAPIDA)")
print("="*70)
print(f"Numero totale frame:  {num_steps}")
print(f"Durata simulazione:   {time[-1]:.2f} s\n")

print("--- VALORI MASSIMI E MINIMI ---")
print("Se Min e Max sono 0.0, il sensore non sta leggendo nulla!")
print("Se il Flap tocca i 30.0, il PID e' andato in saturazione!")
print(f"- Rollio Base (990000):   Min = {np.min(base_E1):>7.2f}° | Max = {np.max(base_E1):>7.2f}°")
print(f"- Rollio Ala (990001):    Min = {np.min(ala_E1):>7.2f}° | Max = {np.max(ala_E1):>7.2f}°")
print(f"- Velocità Rollio:        Min = {np.min(base_omega1):>7.2f}°/s| Max = {np.max(base_omega1):>7.2f}°/s")
print(f"- Comando PID (Flap):     Min = {np.min(flap_cmd_deg):>7.2f}° | Max = {np.max(flap_cmd_deg):>7.2f}°\n")

print("--- ANALISI PASSO-PASSO (Primi 5 frame + Ultimo) ---")
print(f"{'Tempo [s]':<10} | {'Roll Base [°]':<14} | {'Roll Ala [°]':<13} | {'Omega [°/s]':<13} | {'PID Flap [°]':<13}")
print("-" * 70)
for i in range(min(5, num_steps)):
    print(f"{time[i]:<10.3f} | {base_E1[i]:<14.3f} | {ala_E1[i]:<13.3f} | {base_omega1[i]:<13.3f} | {flap_cmd_deg[i]:<13.3f}")

print("   ...     |      ...       |      ...      |      ...      |      ...")
print(f"{time[-1]:<10.3f} | {base_E1[-1]:<14.3f} | {ala_E1[-1]:<13.3f} | {base_omega1[-1]:<13.3f} | {flap_cmd_deg[-1]:<13.3f}")
print("="*70 + "\n")

# ==============================================================================
# PLOTTING DEI DATI: DASHBOARD COMPLETA
# ==============================================================================
print("Generazione dashboard grafica...")
plt.style.use('seaborn-v0_8-darkgrid')
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

ax1.plot(time, base_E1, 'b-', linewidth=2.5, label=f'Rollio Baricentro ({ID_BASE})')
ax1.plot(time, ala_E1, 'c--', linewidth=2, label=f'Rollio Ala ({ID_ALA})')
ax1.axhline(SETPOINT_ROLLIO, color='r', linestyle='--', linewidth=2, label='Target (10°)')
ax1.fill_between(time, base_E1, SETPOINT_ROLLIO, color='red', alpha=0.1, label='Errore Proporzionale')
ax1.set_ylabel('Rollio [Gradi]', fontsize=11, fontweight='bold')
ax1.set_title('Dinamica di Rollio vs Flessione Alare', fontsize=14, fontweight='bold')
ax1.legend(loc='lower right')

ax2.plot(time, base_E2, 'r-', linewidth=2, label='Beccheggio (E2)')
ax2.plot(time, base_E3, 'g-', linewidth=2, label='Imbardata (E3)')
ax2.plot(time, ala_E2, 'orange', linestyle='--', label='Beccheggio Ala')
ax2.set_ylabel('Gradi', fontsize=11, fontweight='bold')
ax2.set_title('Effetti Secondari (Moti Accoppiati)', fontsize=12, fontweight='bold')
ax2.legend(loc='upper right')

ax3.plot(time, base_omega1, 'm-', linewidth=2.5, label='Rateo di Rollio (Omega 1)')
ax3.axhline(0.0, color='k', linestyle='-', linewidth=1)
ax3.set_ylabel('Velocità [Gradi/s]', fontsize=11, fontweight='bold')
ax3.set_title('Velocità Angolare', fontsize=12, fontweight='bold')
ax3.legend(loc='upper right')

ax4.plot(time, flap_cmd_deg, 'g-', linewidth=2.5, label='Comando Flap (Output PID)')
ax4.axhline(30.0, color='r', linestyle=':', linewidth=2, label='Limite Superiore (+30°)')
ax4.axhline(-30.0, color='r', linestyle=':', linewidth=2, label='Limite Inferiore (-30°)')
ax4.set_xlabel('Tempo [s]', fontsize=12, fontweight='bold')
ax4.set_ylabel('Deflessione [Gradi]', fontsize=11, fontweight='bold')
ax4.set_title('Sforzo del Controller PID', fontsize=12, fontweight='bold')
ax4.legend(loc='lower right')

plt.tight_layout()

nome_output = 'dashboard_aeroelastica_pid.png'
plt.savefig(nome_output, dpi=150)
print(f"Finito! Immagine '{nome_output}' salvata.")

try:
    plt.show()
except Exception:
    pass