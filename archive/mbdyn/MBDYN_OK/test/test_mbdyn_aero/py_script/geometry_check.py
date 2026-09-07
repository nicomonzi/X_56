import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# =========================================================================
# VISUALIZZATORE GEOMETRIA X-56A - CORREZIONE PARALLELISMO ALA/TRAVE
# =========================================================================

# --- 1. DATI GEOMETRICI (da dustmesh.txt e Nastran) ---
# Dati estratti dalle prime 5 regioni (fino a Y=50, dove la trave è dritta)
Y_breaks = [0.0, 3.038, 9.350, 19.021, 25.665, 50.000]
C_breaks = [90.53, 86.06, 72.49, 52.25, 43.88, 23.94]
X_LE_breaks = [0.0, 3.3176, 14.2369, 30.5949, 36.2774, 46.6224]

# Coordinate esatte dei nodi estratti dal file Nastran
nodes_y = [
    0.0, 11.04, 22.07, 31.38, 40.69, 50.00, 58.43, 66.86, 75.29, 83.71, 
    92.14, 100.57, 109.00, 117.43, 125.86, 134.29, 142.71, 151.14, 159.57, 168.00
]
nodes_x_nastran = [
    160.15, 160.15, 160.15, 160.15, 160.15, 160.15, 163.46, 166.78, 170.10, 173.42, 
    176.74, 180.06, 183.38, 186.70, 190.02, 193.33, 196.65, 199.97, 203.29, 206.61
]

# Shift per allineare l'apice del DUST (LE=0) con il sistema Nastran
X_SHIFT = 160.15 - 58.38 

# --- 2. EQUAZIONE PROFILO NACA 0012 ---
def naca0012_half_thickness(x_c, chord):
    t = 0.12 # Spessore 12%
    yt = 5 * t * chord * (0.2969 * np.sqrt(x_c) - 0.1260 * x_c - 
                          0.3516 * (x_c**2) + 0.2843 * (x_c**3) - 0.1015 * (x_c**4))
    return yt

# --- 3. FUNZIONE DI INTERPOLAZIONE (Corretta per parallelismo) ---
def get_section(y):
    y = abs(y)
    
    # 1. Posizione esatta della trave secondo il Nastran
    beam_nastran = np.interp(y, nodes_y, nodes_x_nastran)
    beam = beam_nastran - X_SHIFT
    
    if y <= 50.0:
        # Sezione centrale rigida (Centerbody): seguiamo la mesh DUST complessa
        le = np.interp(y, Y_breaks, X_LE_breaks)
        c = np.interp(y, Y_breaks, C_breaks)
        te = le + c
    else:
        # Ala esterna flessibile (Y > 50): FORZIAMO IL PARALLELISMO CON LA TRAVE
        # Calcoliamo l'offset trave-LE esattamente al nodo di giunzione (Y=50)
        beam_at_50 = 160.15 - X_SHIFT
        le_at_50 = 46.6224
        dist_beam_le = beam_at_50 - le_at_50
        
        # Manteniamo questa distanza costante verso l'esterno
        le = beam - dist_beam_le
        c = 23.94 # La corda è costante a 23.94 nell'ala esterna
        te = le + c
        
    return le, te, beam

# --- 4. GENERAZIONE DATI PER IL PLOT ---
X_LE, Y_LE, Z_LE = [], [], []
X_TE, Y_TE, Z_TE = [], [], []
X_Beam, Y_Beam, Z_Beam = [], [], []
X_AC, Y_AC, Z_AC = [], [], []

all_y = [-y for y in reversed(nodes_y)] + nodes_y[1:]

# --- 5. PLOTTING ---
fig = plt.figure(figsize=(14, 9))
ax = fig.add_subplot(111, projection='3d')
s = np.linspace(0, 1, 30)

for y in all_y:
    le, te, beam = get_section(y)
    chord = te - le
    ac = le + 0.25 * chord
    
    Y_LE.append(y); X_LE.append(le); Z_LE.append(0)
    Y_TE.append(y); X_TE.append(te); Z_TE.append(0)
    Y_Beam.append(y); X_Beam.append(beam); Z_Beam.append(0)
    Y_AC.append(y); X_AC.append(ac); Z_AC.append(0)
    
    x_profile = le + s * chord
    z_upper = naca0012_half_thickness(s, chord)
    z_lower = -z_upper
    y_profile = np.full_like(x_profile, y)
    
    ax.plot(x_profile, y_profile, z_upper, color='#555555', alpha=0.5, linewidth=1)
    ax.plot(x_profile, y_profile, z_lower, color='#555555', alpha=0.5, linewidth=1)

# Plot Linee guida
ax.plot(X_LE, Y_LE, Z_LE, 'r-', linewidth=2.5, label='Leading Edge (Parallelo)')
ax.plot(X_TE, Y_TE, Z_TE, 'g-', linewidth=2.5, label='Trailing Edge (Parallelo)')
ax.plot(X_Beam, Y_Beam, Z_Beam, 'b-o', linewidth=2, markersize=4, label='Asse Elastico (Nastran)')
ax.plot(X_AC, Y_AC, Z_AC, 'y--', linewidth=2, label='Centro Aerodinamico (25% c)')

# Impostazioni visuali
ax.set_title("Geometria X-56A: Correzione Parallelismo Trave-Ala", fontsize=14, fontweight='bold', pad=20)
ax.set_xlabel('Asse X (Corda) [in]')
ax.set_ylabel('Asse Y (Apertura) [in]')
ax.set_zlabel('Asse Z (Spessore) [in]')
ax.invert_xaxis() 

x_min, x_max = min(X_LE), max(X_TE)
y_min, y_max = min(Y_LE), max(Y_LE)
z_max = max(C_breaks) * 0.12
z_min = -z_max

ax.set_xlim([x_min, x_max])
ax.set_ylim([y_min, y_max])
ax.set_zlim([z_min, z_max])
ax.set_box_aspect([x_max - x_min, y_max - y_min, (z_max - z_min) * 3.0]) 

ax.view_init(elev=25, azim=-45)
ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))

plt.tight_layout()
plt.show()