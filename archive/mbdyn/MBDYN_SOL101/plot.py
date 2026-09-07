import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

# ==========================================================
# IMPOSTAZIONI
# ==========================================================
file_nc = '/home/nicomonzi/TESI/SOL101/MBDYN_SOL101/xxx.nc' 

# Inserisci qui gli ID dei nodi dove hai applicato il carico
nodi_interesse = ['990020', '991020'] 
# ==========================================================

try:
    data = nc.Dataset(file_nc, 'r')
except FileNotFoundError:
    print(f"Errore: Il file '{file_nc}' non è stato trovato.")
    exit()

# Estrae il vettore del tempo
try:
    time = data.variables['time'][:]
except KeyError:
    print("Errore: Variabile 'time' non trovata nel file NetCDF.")
    data.close()
    exit()

# Prepara la figura
plt.figure(figsize=(10, 6))

colori = ['b', 'r', 'g', 'm', 'c']


node_1 = 'node.struct.990020.X'
node_2 = 'node.struct.991020.X'
    
# Estrae le posizioni assolute [X, Y, Z] per tutti gli istanti di tempo
pos_1 = data.variables[node_1][:]
pos_2 = data.variables[node_2][:]
        
delta_z1 = pos_1[:, 2] - pos_1[0, 2] 
delta_z2 = pos_2[:, 2] - pos_2[0, 2] 

# Plotta lo spostamento lungo l'asse Z
plt.plot(time, delta_z1, label='Nodo 990020 (Tip)', linewidth=2)
plt.plot(time, -delta_z2, label='Nodo 991020 (Tip)', linewidth=2)

data.close()

# ==========================================================
# PERSONALIZZAZIONE GRAFICO
# ==========================================================
plt.title('Risposta Dinamica: Spostamento Verticale (Asse Z) ai Nodi Caricati', fontsize=14)
plt.xlabel('Tempo (s)', fontsize=12)
plt.ylabel('Deformazione Z (in)', fontsize=12)

# Aggiunge una griglia per facilitare la lettura dei valori
plt.grid(True, linestyle='--', alpha=0.7)

# Evidenzia lo zero
plt.axhline(0, color='black', linewidth=1)

plt.legend(fontsize=12)
plt.tight_layout()

# Mostra il plot a schermo
plt.show()