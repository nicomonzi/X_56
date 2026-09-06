import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from netCDF4 import Dataset
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==========================================
# 1. PARAMETRI E DATI DEL MODELLO
# ==========================================
NETCDF_FILE = '/home/nicomonzi/TESI/test/test_mbdyn_aero/output/x.nc' 

# Nodi STRUTTURALI completi (inclusi quelli senza aerodinamica 21, 22, 23)
struct_nodes_left = list(range(990001, 990024))   # Da 990001 a 990023
struct_nodes_right = list(range(991002, 991024))  # Da 991002 a 991023
all_nodes = struct_nodes_left + struct_nodes_right

# Dati AERODINAMICI (si fermano al nodo 20)
aero_data_left = {
    990001: (11.0400, 90.5300, -35.7475),
    990002: (11.0350, 68.9531, -24.0462),
    990003: (10.1700, 48.4089, -13.0752),
    990004: (9.3100, 39.1972, -9.8739),
    990005: (9.3100, 31.5686, -7.8232),
    990006: (8.8700, 23.9400, -5.7726),
    **{i: (8.4300, 23.9400, -5.7726) for i in range(990007, 990020)},
    990020: (4.2150, 23.9400, -5.7726)
}

aero_data_right = {
    991002: (11.0350, 68.9531, -24.0462),
    991003: (10.1700, 48.4089, -13.0752),
    991004: (9.3100, 39.1972, -9.8739),
    991005: (9.3100, 31.5686, -7.8232),
    991006: (8.8700, 23.9400, -5.7726),
    **{i: (8.4300, 23.9400, -5.7726) for i in range(991007, 991020)},
    991020: (4.2150, 23.9400, -5.7726)
}

# ==========================================
# 2. FUNZIONI DI LETTURA E GEOMETRIA
# ==========================================
def load_mbdyn_kinematics(filepath, nodes):
    try:
        nc = Dataset(filepath, 'r')
        time = nc.variables['time'][:]
        pos_dict = {}
        for node in nodes:
            var_name = f'node.struct.{node}.X'
            if var_name in nc.variables:
                pos_dict[node] = nc.variables[var_name][:]
            else:
                print(f"Attenzione: Variabile {var_name} non trovata.")
        nc.close()
        return time, pos_dict
    except Exception as e:
        print(f"Errore NetCDF: {e}")
        return None, None

def get_naca0012_faces(pos, aero_props, n_points=10):
    """
    Genera le facce 3D di un profilo NACA0012 per l'elemento aerodinamico.
    """
    x, y, z = pos
    width, chord, offset = aero_props
    
    x_ac = x + offset
    x_le = x_ac - 0.25 * chord  # Bordo d'attacco
    
    y_in = y - width / 2
    y_out = y + width / 2
    
    # Creiamo i punti lungo la corda (da 0 a 1). 
    # Usiamo il coseno per addensare i punti vicino a LE e TE per una curva più morbida.
    beta = np.linspace(0, np.pi, n_points)
    x_c = 0.5 * (1 - np.cos(beta)) 
    
    # Formula esatta per lo spessore del profilo NACA0012
    t = 0.12
    yt = 5 * t * (0.2969 * np.sqrt(x_c) - 0.1260 * x_c - 0.3516 * (x_c**2) + 0.2843 * (x_c**3) - 0.1015 * (x_c**4))
    
    faces = []
    
    # Costruiamo il guscio superiore e inferiore a strisce
    for i in range(n_points - 1):
        x1 = x_le + x_c[i] * chord
        x2 = x_le + x_c[i+1] * chord
        
        # Dorso (Upper surface)
        z1_up = z + yt[i] * chord
        z2_up = z + yt[i+1] * chord
        faces.append([
            [x1, y_in, z1_up], [x2, y_in, z2_up], [x2, y_out, z2_up], [x1, y_out, z1_up]
        ])
        
        # Ventre (Lower surface)
        z1_lo = z - yt[i] * chord
        z2_lo = z - yt[i+1] * chord
        faces.append([
            [x1, y_in, z1_lo], [x2, y_in, z2_lo], [x2, y_out, z2_lo], [x1, y_out, z1_lo]
        ])
        
    return faces

# ==========================================
# 3. SETUP DELL'ANIMAZIONE E LIMITI DINAMICI
# ==========================================
time, positions = load_mbdyn_kinematics(NETCDF_FILE, all_nodes)

if time is None or not positions:
    print("Esco per errore di lettura dati.")
    exit()

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Linee strutturali modificate (ora includono tutti i nodi)
line_left, = ax.plot([], [], [], 'ko-', lw=2, markersize=5, label='Trave Sinistra (Completa)')
line_right, = ax.plot([], [], [], 'bo-', lw=2, markersize=5, label='Trave Destra (Completa)')

# Pannelli aerodinamici (colore azzurro semitrasparente, bordi sottili)
panels_collection = Poly3DCollection([], alpha=0.4, facecolors='cyan', edgecolors='blue', linewidths=0.2)
ax.add_collection3d(panels_collection)

# Calcolo dinamico dei limiti
all_x, all_y, all_z = [], [], []
for n in positions:
    all_x.extend(positions[n][:, 0])
    all_y.extend(positions[n][:, 1])
    all_z.extend(positions[n][:, 2])

min_x, max_x = np.min(all_x), np.max(all_x)
min_y, max_y = np.min(all_y), np.max(all_y)
min_z, max_z = np.min(all_z), np.max(all_z)

lim_x = (min_x - 10, max_x + 10)
lim_y = (min_y - 10, max_y + 10)
lim_z = (min_z - 30, max_z + 30)

ax.set_xlim(lim_x)
ax.set_ylim(lim_y)
ax.set_zlim(lim_z)

# Proporzioni 1:1:1
ax.set_box_aspect([lim_x[1]-lim_x[0], lim_y[1]-lim_y[0], lim_z[1]-lim_z[0]])

ax.set_xlabel('X (Streamwise) [in]')
ax.set_ylabel('Y (Spanwise) [in]')
ax.set_zlabel('Z (Vertical) [in]')
ax.legend()

def update(frame):
    # Aggiorna linea sinistra usando la lista completa dei nodi STRUTTURALI
    x_l = [positions[n][frame, 0] for n in struct_nodes_left if n in positions]
    y_l = [positions[n][frame, 1] for n in struct_nodes_left if n in positions]
    z_l = [positions[n][frame, 2] for n in struct_nodes_left if n in positions]
    line_left.set_data(x_l, y_l)
    line_left.set_3d_properties(z_l)
    
    # Aggiorna linea destra
    x_r = [positions[n][frame, 0] for n in struct_nodes_right if n in positions]
    y_r = [positions[n][frame, 1] for n in struct_nodes_right if n in positions]
    z_r = [positions[n][frame, 2] for n in struct_nodes_right if n in positions]
    line_right.set_data(x_r, y_r)
    line_right.set_3d_properties(z_r)

    # Aggiorna il volume aerodinamico NACA0012
    verts = []
    for n in aero_data_left.keys():
        if n in positions:
            # Usiamo extend invece di append perché get_naca0012_faces restituisce una LISITA di facce
            verts.extend(get_naca0012_faces(positions[n][frame], aero_data_left[n]))
            
    for n in aero_data_right.keys():
        if n in positions:
            verts.extend(get_naca0012_faces(positions[n][frame], aero_data_right[n]))

    panels_collection.set_verts(verts)
    
    ax.set_title(f'Animazione MBDyn (NACA0012) - Tempo: {time[frame]:.3f} s')
    
    return line_left, line_right, panels_collection

# Creazione dell'animazione
ani = animation.FuncAnimation(fig, update, frames=len(time), interval=20, blit=False)

plt.show()