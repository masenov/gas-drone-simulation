import numpy as np
import GPy


from IPython.display import display
GPy.plotting.change_plotting_library('matplotlib')

X = np.random.uniform(-3.,3.,(20,1))
Y = np.sin(X) + np.random.randn(20,1)*0.05

kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)

m = GPy.models.GPRegression(X,Y,kernel)


display(m)
fig = m.plot(plot_density=True)
fig.figure.savefig("test.pdf")


m.optimize(messages=True)
m.optimize_restarts(num_restarts = 10)

display(m)
fig = m.plot(plot_density=True)
fig.figure.savefig("optimized.pdf")


# sample inputs and outputs
X = np.random.uniform(-3.,3.,(50,2))
Y = np.sin(X[:,0:1]) * np.sin(X[:,1:2])+np.random.randn(50,1)*0.05

# define kernel
ker = GPy.kern.Matern52(2,ARD=True) + GPy.kern.White(2)

# create simple GP model
m = GPy.models.GPRegression(X,Y,ker)

# optimize and plot
m.optimize(messages=True,max_f_eval = 1000)
fig = m.plot()
fig.figure.savefig("2d_example.pdf")
display(m)

input("Press Enter to continue...")
