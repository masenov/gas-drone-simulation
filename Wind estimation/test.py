import numpy as np
import matplotlib.pyplot as plt
import GPy
from IPython.display import display
import plotly 
plotly.tools.set_credentials_file(username='masenov', api_key='xtuT0tL0GZYldyT3JXGU')
GPy.plotting.change_plotting_library('matplotlib')

"""
V = np.array([[1,1],[-2,2],[4,-7]])
direction = np.pi/4
magnitude = 1

for i in range(-5,5):
    for j in range(-5,5):
        plt.quiver([i],[j], np.cos(direction)*magnitude, np.sin(direction)*magnitude, color=['r'], scale=21)


plt.show()

plt.figure()
Xs = np.array([[-2,4,3]])
Ys = np.array([[-1,1,2]])
direction = np.array([np.pi/4, np.pi/3, np.pi/4])
magnitude = np.array([1,2.1,0.5])
plt.quiver(Xs,Ys, np.cos(direction)*magnitude, np.sin(direction)*magnitude, color='r', scale=21)
plt.show()

Xs = Xs.reshape((3,1))
Ys = Ys.reshape((3,1))
X = np.concatenate((Xs,Ys), axis=1)
print (X.shape,direction.shape,magnitude.shape)
"""
# sample inputs and outputs
n = 30
X = np.random.uniform(-5.,5.,(n,2))
Y = np.sin(X[:,0:1]) * np.sin(X[:,1:2]) * 3.14 +np.random.randn(n,1)*0.05
Z = np.cos(X[:,0:1]) * np.cos(X[:,1:2])+np.random.randn(n,1)*0.05
plt.quiver(X[:,0],X[:,1], np.cos(Y)*Z, np.sin(Y)*Z, color=['r'], scale=21, width=0.007, headwidth=2)
plt.show()
print (Y)
#Y = direction.reshape((3,1))
#Z = magnitude.reshape((3,1))
# define kernel
ker = GPy.kern.Matern52(2,ARD=True) + GPy.kern.White(2)
ker= GPy.kern.RBF(input_dim=2, variance=2., lengthscale=1.)+ GPy.kern.White(2)

# create simple GP model
m = GPy.models.GPRegression(X,Y,ker)
ker= GPy.kern.RBF(input_dim=2, variance=2., lengthscale=1.)+ GPy.kern.White(2)
n = GPy.models.GPRegression(X,Z,ker)

# optimize and plot
m.optimize(messages=True,max_f_eval = 1000)
n.optimize(messages=True,max_f_eval = 1000)
"""
fig = m.plot()
fig.figure.savefig("2d_example.pdf")
fig = n.plot()
fig.figure.savefig("2d1_example.pdf")
"""
magnitude = 1
for i in np.arange(-5.0,5.0,0.5):
    for j in np.arange(-5.0,5.0,0.5):
        direction = m.predict(np.array([[i,j]]))[0]
        magnitude = n.predict(np.array([[i,j]]))[0]*2
        #print (direction,magnitude,i,j)
        plt.quiver([i],[j], np.cos(direction)*magnitude, np.sin(direction)*magnitude, color=['r'], scale=21, width=0.007, headwidth=2)
plt.show()
