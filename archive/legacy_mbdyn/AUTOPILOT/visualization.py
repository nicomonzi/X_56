import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks

# --- 1. CONFIGURAZIONE E LETTURA DATI ---
FILE_NC = "/home/nicomonzi/TESI/AUTOPILOT/clamp.nc"
ds = nc.Dataset(FILE_NC, 'r')

t = ds.variables['time'][:]
dt = t[1] - t[0]
N = len(t)

# Geometria della trave (ala completa)
nodi_sx = [990023 - i for i in range(23)]
nodi_dx = [991002 + i for i in range(22)]
ordine_nodi = nodi_sx + nodi_dx

lista_posizioni = []
for nodo in ordine_nodi:
    nome_var = f'node.struct.{nodo}.X'
    if nome_var in ds.variables:
        lista_posizioni.append(ds.variables[nome_var][:])
dati_trave = np.stack(lista_posizioni, axis=1)

# Estrazione Spostamento Z (Flessione) ai due Tip
z_sx = ds.variables['node.struct.990020.X'][:, 2]
z_dx = ds.variables['node.struct.991020.X'][:, 2]

# Estrazione Rotazione Y (Torsione) - Gestione dinamica dell'output MBDyn (.Phi o .E)
if 'node.struct.990023.Phi' in ds.variables:
    rot_sx = ds.variables['node.struct.990023.Phi'][:, 1]  # Componente Y del vettore di rotazione
    rot_dx = ds.variables['node.struct.991023.Phi'][:, 1]
elif 'node.struct.990023.E' in ds.variables:
    rot_sx = ds.variables['node.struct.990023.E'][:, 1]    # Componente Y in angoli di Eulero
    rot_dx = ds.variables['node.struct.991023.E'][:, 1]
else:
    print("Mancano le variabili di orientazione (.Phi/.E). Torsione impostata a zero.")
    rot_sx, rot_dx = np.zeros_like(z_sx), np.zeros_like(z_dx)

ds.close()

z_sx_c = z_sx - np.mean(z_sx)
z_dx_c = z_dx - np.mean(z_dx)
theta_sx_c = np.degrees(rot_sx) - np.mean(np.degrees(rot_sx))
theta_dx_c = np.degrees(rot_dx) - np.mean(np.degrees(rot_dx))

xf = fftfreq(N, dt)[:N//2]

yf_z = fft(z_sx_c)
amp_z = 2.0 / N * np.abs(yf_z[:N//2])
freq_dom_z = xf[np.argmax(amp_z)]

yf_theta = fft(theta_sx_c)
amp_theta = 2.0 / N * np.abs(yf_theta[:N//2])
freq_dom_theta = xf[np.argmax(amp_theta)]

picchi_idx, _ = find_peaks(z_sx_c)
smorzamento_zeta = 0.0
if len(picchi_idx) > 1:
    z_picchi = z_sx_c[picchi_idx]
    delta = (1 / (len(z_picchi) - 1)) * np.log(z_picchi[0] / z_picchi[-1])
    smorzamento_zeta = delta / np.sqrt(4 * np.pi**2 + delta**2)

# --- 2. STAMPA INFO NEL TERMINALE ---
coalescenza = abs(freq_dom_z - freq_dom_theta)

testo_output = (
    f"ANALISI DINAMICA STRUTTURALE\n"
    f"• Freq Flessione (Z): {freq_dom_z:.2f} Hz\n"
    f"• Freq Torsione (\u03B8_y): {freq_dom_theta:.2f} Hz\n"
    f"• Smorzamento Flex (\u03B6): {smorzamento_zeta*100:.2f} %\n"
    f"============================================\n"
)
print(testo_output)


# 1
fig1 = plt.figure(figsize=(8, 6))
ax3d = fig1.add_subplot(111, projection='3d')
ax3d.plot(dati_trave[0, :, 0], dati_trave[0, :, 1], dati_trave[0, :, 2], color='lightgray', marker='o', markersize=2, label='Indeformata')
linea_trave, = ax3d.plot([], [], [], 'b-o', linewidth=2, markersize=3, label='Deformata')
ax3d.set_xlabel('X (Corda) [in]')
ax3d.set_ylabel('Y (Apertura) [in]')
ax3d.set_zlabel('Z (Deformata) [in]')
ax3d.invert_xaxis()

X, Y, Z = dati_trave[:,:,0], dati_trave[:,:,1], dati_trave[:,:,2]
max_range = max(X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()) / 2.0
mid_x, mid_y, mid_z = (X.max()+X.min())*0.5, (Y.max()+Y.min())*0.5, (Z.max()+Z.min())*0.5
ax3d.set_xlim(mid_x - max_range, mid_x + max_range)
ax3d.set_ylim(mid_y - max_range, mid_y + max_range)
ax3d.set_zlim(mid_z - max_range, mid_z + max_range)
ax3d.legend(loc='upper left')
fig1.tight_layout()

# 2
fig2 = plt.figure(figsize=(10, 5))
ax_tz = fig2.add_subplot(111)
ax_t_theta = ax_tz.twinx()

p1, = ax_tz.plot(t, z_sx_c, color='tab:blue', label='Flessione Z (Tip Sx)')
p2, = ax_t_theta.plot(t, theta_sx_c, color='tab:red', linestyle='--', label='Torsione $\\theta_y$ (Tip Sx)')

ax_tz.set_title('Risposta Temporale al Tip Sinistro')
ax_tz.set_xlabel('Tempo [s]')
ax_tz.set_ylabel('Spostamento Z [in]', color='tab:blue')
ax_t_theta.set_ylabel('Angolo di Torsione [deg]', color='tab:red')
ax_tz.tick_params(axis='y', labelcolor='tab:blue')
ax_t_theta.tick_params(axis='y', labelcolor='tab:red')

linea_cursore, = ax_tz.plot([t[0], t[0]], [np.min(z_sx_c), np.max(z_sx_c)], 'g:', linewidth=2, zorder=5)
ax_tz.grid(True)
fig2.tight_layout()

fig3 = plt.figure(figsize=(10, 5))
ax_fz = fig3.add_subplot(111)
ax_f_theta = ax_fz.twinx()

f1, = ax_fz.plot(xf, amp_z, color='tab:blue', alpha=0.8, label='FFT Z')
f2, = ax_f_theta.plot(xf, amp_theta, color='tab:red', alpha=0.6, linestyle='--', label='FFT $\\theta_y$')

ax_fz.set_title('Spettri di Frequenza Sovrapposti')
ax_fz.set_xlabel('Frequenza [Hz]')
ax_fz.set_ylabel('Ampiezza Flessione', color='tab:blue')
ax_f_theta.set_ylabel('Ampiezza Torsione', color='tab:red')
ax_fz.set_xlim(0, 20)
ax_fz.grid(True)
fig3.tight_layout()

# 4
SALTA_FRAME = max(1, N // 200)

def update(frame_idx):
    v_idx = frame_idx * SALTA_FRAME
    if v_idx >= N: v_idx = N - 1
    
    # Aggiorna la trave 3D (Figura 1)
    linea_trave.set_data(dati_trave[v_idx, :, 0], dati_trave[v_idx, :, 1])
    linea_trave.set_3d_properties(dati_trave[v_idx, :, 2])
    ax3d.set_title(f"Deformata Trave 3D | t: {t[v_idx]:.3f} s", fontweight='bold')
    
    # Aggiorna la linea verticale del tempo corrente (Figura 2)
    linea_cursore.set_xdata([t[v_idx], t[v_idx]])
    fig2.canvas.draw_idle()  # Forza l'aggiornamento della Figura 2
    
    return linea_trave, linea_cursore

# Alleghiamo l'animazione alla prima figura
ani = animation.FuncAnimation(fig1, update, frames=N//SALTA_FRAME, interval=40, blit=False)

plt.show()