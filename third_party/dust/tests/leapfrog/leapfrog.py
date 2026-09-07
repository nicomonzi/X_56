import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
import concurrent.futures as cf # Import the concurrent.futures library (for parallel processing)
from itertools import repeat    # Import the repeat function from the itertools library 
import h5py 
import pandas as pd 

plt.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 14 

def calc_rings_weighted(leapfrog_file, nrings, intervals):
    
    # Iterate over each ring
    outZ = np.zeros((nrings, 3))
    outR = np.zeros(nrings) 
    for ri in range(nrings):

        # Calculate centroid
        outZ[ri] = np.zeros(3)
        magGammatot = 0
        for pi in range(int(intervals[ri][0]), int(intervals[ri][1])):

            vorticity = leapfrog_file.vorticity[pi, :]
            position = leapfrog_file.position[pi, :] 
            normGamma = np.linalg.norm(vorticity)
            magGammatot += normGamma

            for i in range(3):
                outZ[ri][i] += normGamma*position[i] 

        outZ[ri] /= magGammatot

        # Calculate ring radius #and cross-section radius
        outR[ri] = 0
        for pi in range(int(intervals[ri][0]), int(intervals[ri][1])):
            
            vorticity = leapfrog_file.vorticity[pi, :]
            position = leapfrog_file.position[pi, :] 
            normGamma = np.linalg.norm(vorticity)

            outR[ri] += normGamma*np.sqrt(  (position[0] - outZ[ri][0])**2 + 
                                            (position[1] - outZ[ri][1])**2 + 
                                            (position[2] - outZ[ri][2])**2)

        outR[ri] /= magGammatot

    return outZ, outR

def process_time(t, nrings):

    # Format the filename using the current time step (t)
    filehdf5 = f'output/ring_re3000_res_{str(time[t]+1).zfill(4)}.h5'
    print(filehdf5)  # Print the filename

    # Open the HDF5 file in read mode
    leapfrog_file = h5py.File(filehdf5, 'r')

    # Define the dataset paths in the HDF5 file
    sol_dsets = ['/ParticleWake/WakePoints', '/ParticleWake/WakeVort'] 

    # Load the position data from the HDF5 file into a NumPy array
    leapfrog_file.position = np.array(leapfrog_file[sol_dsets[0]])

    # Load the vorticity data from the HDF5 file into a NumPy array
    leapfrog_file.vorticity = np.array(leapfrog_file[sol_dsets[1]])

    # Get the number of particles from the shape of the position data
    nparticles = leapfrog_file.position.shape[0]

    # Initialize a 2x2 array for the intervals
    intervals = np.zeros((2, 2)) 

    # Set the first interval to [0, nparticles/2] (rounded down)
    intervals[0] = np.array([0, int(np.floor(nparticles/2))])

    # Set the second interval to [nparticles/2 (rounded up), nparticles]
    intervals[1] = np.array([int(np.ceil(nparticles/2)), nparticles])  

    # Call the calc_rings_weighted function with the loaded data and return its result
    return calc_rings_weighted(leapfrog_file, nrings, intervals) 

if __name__ == "__main__":

    time = np.arange(0, 1999)
    outZ = np.zeros((len(time), 2, 3))
    outR = np.zeros((len(time), 2)) 
    nrings = 2 # ring numbers 
    # reference values from 'Leapfrogging of multiple coaxial viscous vortex rings, 2015'  
    ring_1_lb_data = pd.read_csv('ring_1_lb.dat', delim_whitespace=True, header=None)
    ring_1_lb = ring_1_lb_data.to_numpy()
    ring_2_lb_data = pd.read_csv('ring_2_lb.dat', delim_whitespace=True, header=None) 
    ring_2_lb = ring_2_lb_data.to_numpy() 
    
    with cf.ProcessPoolExecutor() as executor: 
        results = executor.map(process_time, range(len(time)), repeat(nrings)) 
        for t, result in enumerate(results): 
            for ri in range(nrings): 
                outR[t, ri] = result[1][ri] 
                outZ[t, ri] = result[0][ri] 
    
    fig = plt.figure(figsize=plt.figaspect(0.5))
    plt.plot(outZ[:, 1, 2]/outR[0,1] + outZ[0, 0, 2], outR[:,1]/outR[0, 1], label='Ring 1', color='red',    linewidth=2.0)
    plt.plot(outZ[:, 0, 2]/outR[0,0] + outZ[0, 0, 2], outR[:,0]/outR[0, 0], label='Ring 2', color='blue',       linewidth=2.0)  
    plt.scatter(ring_1_lb[:,1], ring_1_lb[:,0], label='LB Ring 1', color='black', s=8.0, marker = 'o') 
    plt.scatter(ring_2_lb[:,1], ring_2_lb[:,0], label='LB Ring 2', color='black', s=8.0, marker = 'd')  
    # add labels
    plt.xlabel(r'Ring centroid $\frac{Z}{R_0}$')
    plt.ylabel(r'Ring radius   $\frac{R}{R_0}$')  
    # Set the limits of the x and y axes
    plt.xlim(0, 8)
    plt.ylim(0.6, 1.4)
    plt.legend()  
    plt.grid(True)
    plt.show()  
    fig.savefig('leapfrog3000.png')
    plt.close(fig) 