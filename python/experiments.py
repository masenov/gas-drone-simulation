from fluidsimulation import *
from visualize import *

def run(fs,f,timesteps):
    for i in range(timesteps):
        fs.update(f)

    Visualize(filename, timesteps)
    plt.close()


timesteps = 100
filename = 'mytestfile'



fs = FluidSimulation(256)
#fs.resetVelocity()

f = h5py.File(filename + ".hdf5", "w")



run(fs,f,timesteps)
