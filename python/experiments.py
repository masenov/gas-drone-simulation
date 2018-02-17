from fluidsimulation import *
from visualize import *
import numpy as np
import time



def sample(f, size=1000, timestep=40, locations=5, gridSize=130, offset=20):
    nx, ny = (locations, locations)
    x = np.linspace(offset, gridSize-offset, nx)
    y = np.linspace(offset, gridSize-offset, ny)
    xv, yv = np.meshgrid(x, y)
    xv = xv.astype(int)
    yv = yv.astype(int)
    print (xv.astype(int))
    # For wind 0/1 for wind speed and direction
    # f.file['wind0'][0][1]
    # f.file['test0'][64][64]
    n = size//timestep
    concetrations = np.zeros((n,locations,locations))
    wind = np.zeros((n,2))
    for i in range(n):
        for j in range(locations):
            for k in range(locations):
                concetrations[i][j][k] = f.file['test'+str(i*40)][yv[j][k]][xv[j][k]]
        for j in range(2):
            wind[i][j] = f.file['wind'+str(i*40)][0][j]
    return concetrations, wind


def sample_vis(concentrations):
    plt.figure()
    plt.ion()     # turns on interactive mode
    for i in range(concentrations.shape[0]):
        plt.imshow(concentrations[i,:,:]/concentrations.max())
        plt.show()
        print concentrations
        plt.pause(0.01)
    plt.ioff()     # turns on interactive mode


def real_sample_vis(concentrations):
    plt.ion()     # turns on interactive mode
    true_sample = np.zeros((5,5))
    n = 0
    for i in range(5):
        for j in range(5):
            true_sample[i,j] = concentrations[n,i,j]
            n = n + 1

    plt.imshow(true_sample)
    plt.show()
    plt.ioff()     # turns on interactive mode




def run(fs,f,timesteps, simulate=True):
    if (simulate):
        for i in range(timesteps):
            fs.update(f)

    Visualize(filename, timesteps)
    plt.close()


timesteps = 1000
filename = 'mytestfile4'
simulate = False


fs = FluidSimulation(128, windDirection = 60.0)

if (simulate):
    f = h5py.File(filename + ".hdf5", "w")
else:
    f = h5py.File(filename + ".hdf5", "r")

concentrations, wind = sample(f)
real_sample_vis(concentrations)
sample_vis(concentrations)
run(fs,f,timesteps,simulate=simulate)
