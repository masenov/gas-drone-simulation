from gaussian import *



sigma_v = 0
l = 4
sigma_f = 10

samples = 1
ndim = 50
filename = 'data/data_' + str(sigma_v) + '_' + str(l) + '_' + str(sigma_f) + '_' + str(ndim)+'.out'

if os.path.exists(filename):
    (x,y,z) = loadData(filename)
else:
    cov,x,y = sigma2(ndim,length=10,uniform=False,sigma_v=sigma_v,l=l,sigma_f=sigma_f)
    mu=np.zeros(ndim*ndim)
    z=np.random.multivariate_normal(mu, cov, samples).T
    saveData(x,y,z,filename)


z = z.reshape(ndim,ndim).T
x = x.reshape(ndim,ndim).T
y = y.reshape(ndim,ndim).T
print ("All")
print (x,y,z)
print ("X")
print (x)
print ("Y")
print (y)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(x.flatten(),y.flatten(),z.flatten(), cmap=cm.jet)
fig.show()

fig = plt.figure()
surf = plt.contour(x,y,z,cmap=cm.jet)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
