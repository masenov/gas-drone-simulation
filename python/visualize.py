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

fig = plt.figure()


file = h5py.File("mytestfile2.hdf5", "r")
def f(x, y):
    return np.sin(x) + np.cos(y)

x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
i = 0
im = plt.imshow(f(x, y), animated=True)


def updatefig(*args):
    global x, y, i
    im.set_array(file['test'+str(i)])
    i += 1
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=True)
plt.show()
