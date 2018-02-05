"""
=================
An animated image
=================

This example demonstrates how to animate an image.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import sys, traceback



class Visualize:
    def __init__(self, filename, timesteps):
        self.file= h5py.File(filename+".hdf5", "r")
        self.i = 0
        self.fig = plt.figure()
        self.im = plt.imshow(np.random.random((100,100)), animated=True)
        ani = animation.FuncAnimation(self.fig, self.updatefig, frames=timesteps - 2, interval=0, blit=True, repeat=False)
        plt.show(block=False)


    def updatefig(self,*args):
        self.im.set_array(self.file['test'+str(self.i)])
        self.i += 1
        return self.im,


