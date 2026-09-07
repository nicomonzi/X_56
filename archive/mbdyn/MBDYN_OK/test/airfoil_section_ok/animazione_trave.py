import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# =========================================================================
# ANIMAZIONE AEROELASTICA (MBDyn .mov) - NASA X-56A
# =========================================================================

# Nome del file generato da MBDyn
NOME_FILE = "/home/nicomonzi/TESI/test/airfoil_section_ok/output/x.mov" 

# Nodi in ordine per disegnare la linea continua (Da Tip Sinistra a Tip Destra)
nodi_sx = [str(990023 - i) for i in range(23)] # Punta SX -> Radice
nodo_radice = '990000' # Opzionale, se lo hai salvato
nodi_dx = [str(991002 + i) for i in range(22)] # Radice -> Punta DX

# Usiamo una lista per mantenere l'ordine logico di disegno
ordine_nodi = nodi_sx + nodi_dx
tutti_i_nodi = set(ordine_nodi)

# Struttura dati: lista di dizionari (un dizionario = un time-step/fotogramma)
fotogrammi = []
frame_corrente = {}

print(f"Lettura e raggruppamento dei time-step dal file {NOME_FILE}...")

try:
    with open(NOME_FILE, "r") as f:
        for riga in f:
            if not riga.strip(): continue
            dati = riga.split()
            nodo_id = dati[0]
            
            if nodo_id in tutti_i_nodi:
                # Se il nodo è già nel frame corrente, significa che è iniziato un nuovo time-step!
                if nodo_id in frame_corrente:
                    fotogrammi.append(frame_corrente)
                    frame_corrente = {} # Svuotiamo per il nuovo frame
                
                frame_corrente[nodo_id] = [float(dati[1]), float(dati[2]), float(dati[3])]
                
        # Aggiungiamo l'ultimo frame rimasto appeso
        if frame_corrente:
            fotogrammi.append(frame_corrente)

    num_frames = len(fotogrammi)
    print(f"Trovati {num_frames} step temporali.")

    if num_frames == 0:
        print("Errore: Nessun dato valido trovato.")
        exit()

    # --- PREPARAZIONE GRAFICO 3D ---
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Oggetto linea che verrà aggiornato ad ogni frame
    linea_ala, = ax.plot([], [], [], 'b-o', linewidth=2.5, markersize=4)
    titolo_tempo = ax.set_title("Inizializzazione...", fontweight='bold')

    # Estraiamo i limiti globali per fissare la "gabbia" 3D e non farla saltellare
    all_x = [pos[0] for f in fotogrammi for pos in f.values()]
    all_y = [pos[1] for f in fotogrammi for pos in f.values()]
    all_z = [pos[2] for f in fotogrammi for pos in f.values()]

    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)
    min_z, max_z = min(all_z), max(all_z)

    # Margini visivi
    ax.set_xlim(min_x - 10, max_x + 10)
    ax.set_ylim(min_y - 10, max_y + 10)
    ax.set_zlim(min_z - 30, max_z + 30) # Più margine su Z per vedere la flessione

    ax.set_xlabel('Asse X (Corda) [in]')
    ax.set_ylabel('Asse Y (Apertura) [in]')
    ax.set_zlabel('Asse Z (Deformata) [in]')
    ax.invert_xaxis() # Bordo d'attacco in alto

    # Proporzioni 1:1:1
    x_range, y_range, z_range = max_x - min_x, max_y - min_y, max_z - min_z
    ax.set_box_aspect([x_range, y_range, max(z_range, y_range*0.2)])

    # --- FUNZIONE DI AGGIORNAMENTO ANIMAZIONE ---
    # Se hai 10.000 step, l'animazione sarà lenta. Saltiamo qualche frame.
    SALTA_FRAME = max(1, num_frames // 400) # Punta a max ~400 frame a video

    def update(frame_idx):
        vero_idx = frame_idx * SALTA_FRAME
        if vero_idx >= num_frames: vero_idx = num_frames - 1
        
        dati_step = fotogrammi[vero_idx]
        
        # Estraiamo le coordinate nell'ordine corretto da punta a punta
        x_attuali = [dati_step[n][0] for n in ordine_nodi if n in dati_step]
        y_attuali = [dati_step[n][1] for n in ordine_nodi if n in dati_step]
        z_attuali = [dati_step[n][2] for n in ordine_nodi if n in dati_step]

        linea_ala.set_data(x_attuali, y_attuali)
        linea_ala.set_3d_properties(z_attuali)
        
        titolo_tempo.set_text(f"Deformazione Dinamica X-56A | Step: {vero_idx}/{num_frames}")
        return linea_ala, titolo_tempo

    # --- LANCIO ANIMAZIONE ---
    frame_totali = num_frames // SALTA_FRAME
    print(f"Creazione animazione a schermo ({frame_totali} fotogrammi)...")
    
    ani = animation.FuncAnimation(fig, update, frames=frame_totali, interval=30, blit=False)
    
    plt.tight_layout()
    plt.show()

except FileNotFoundError:
    print(f"Errore: Impossibile trovare il file {NOME_FILE}.")