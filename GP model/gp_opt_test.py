#import GPy
import GPyOpt
from numpy.random import seed

def myf(x):
    return (2*x)**2

bounds = [{'name': 'var_1', 'type': 'continuous', 'domain': (-1,1)}]

max_iter = 15

myProblem = GPyOpt.methods.BayesianOptimization(myf,bounds)

myProblem.run_optimization(max_iter)
print (myProblem.x_opt)
print (myProblem.fx_opt)
