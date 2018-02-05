from fluidsimulation import *
from visualize import *


timesteps = 100
filename = 'mytestfile'



fs = FluidSimulation(256)
fs.resetVelocity()

f = h5py.File(filename + ".hdf5", "w")




for i in range(timesteps):
    fs.update(f)

Visualize(filename, timesteps)
plt.close()


for i in range(timesteps):
    fs.update(f)

Visualize(filename, timesteps)
plt.close()
