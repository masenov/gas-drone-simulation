from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import holoviews as hv
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import os.path

interval = 10
interval_size = 0.01

mpl.rcParams.update({'font.size':14})
mpl.rcParams['xtick.major.pad']='8'
mpl.rcParams['ytick.major.pad']='8'
mpl.rcParams['lines.linewidth'] = 4
plt.rc('text', usetex=True)


def gaussian(x,mu,sigma):
    y = np.exp(-(x-mu)**2/(2*sigma**2))/(2*pi*sigma**2)**(1.0/2)
    return y


def pdf_multivariate_gauss(x, mu, cov):
    '''
    Caculate the multivariate normal density (pdf)

    Keyword arguments:
        x = numpy array of a "d x 1" sample vector
        mu = numpy array of a "d x 1" mean vector
        cov = "numpy array of a d x d" covariance matrix
    '''
    assert(mu.shape[0] > mu.shape[1]), 'mu must be a row vector'
    assert(x.shape[0] > x.shape[1]), 'x must be a row vector'
    assert(cov.shape[0] == cov.shape[1]), 'covariance matrix must be square'
    assert(mu.shape[0] == cov.shape[0]), 'cov_mat and mu_vec must have the same dimensions'
    assert(mu.shape[0] == x.shape[0]), 'mu and x must have the same dimensions'
    part1 = 1 / ( ((2* np.pi)**(len(mu)/2)) * (np.linalg.det(cov)**(1/2)) )
    part2 = (-1/2) * ((x-mu).T.dot(np.linalg.inv(cov))).dot((x-mu))
    return float(part1 * np.exp(part2))

def two_dimen_gaussian(mu,cov,x = np.arange(-7, 7, 0.1),y = np.arange(-7, 7, 0.1), samples=False):
    fig_size = 5
    #plt.figure(figsize=[fig_size+1.5, fig_size])
    plt.figure(figsize=[fig_size, fig_size])
    xx, yy = np.meshgrid(x, y, sparse=True)
    z = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            z[i,j] = pdf_multivariate_gauss(np.array([[xx[0][i]], [yy[j][0]]]), mu, cov)
    plt.contourf(x,y,z,cmap=cm.gray)
    x, y = np.random.multivariate_normal([mu[0][0],mu[1][0]], cov, 100).T
    if samples:
        plt.plot(x, y, 'x')
        plt.axis('equal')
        plt.xlim([-4,4])
        plt.ylim([-4,4])





def k(i,j,sigma_f,l):
    return np.power(sigma_f,2)*np.exp(-np.power((i-j),2)/(2*np.power(l,2)))

def direct_delta(i,j):
    if (i==j):
        return 1
    else:
        return 0

def sigma(ndim=6,length=6,uniform=False,sigma_v=0.1,l=3,sigma_f=1):
    if (uniform):
        x = np.random.uniform(1,length,ndim)
        x.sort()
    else:
        x = np.linspace(1, length, ndim)
    cov = np.zeros((ndim,ndim))
    for i in range(ndim):
        for j in range(ndim):
            cov[i,j]=k(x[i],x[j],sigma_f,l)+np.power(sigma_v,2)*direct_delta(i,j)
    return cov,x





def plotGP(x=None,y=None,dots=False,usecolors=True,lw=3,filename=None,xname="X",yname="Y",xlim=None,ylim=None,data=None):
    if (data is not None):
        plt.figure()
        plt.plot(data[0],data[1],'*',c='b',markersize=10)
        plt.xlabel(xname)
        plt.ylabel(yname)
        if (xlim):
            plt.xlim(0.5,0.5+xlim)
        if (ylim):
            plt.ylim(-ylim,ylim)
        plt.savefig(str(filename)+'_data.png', dpi=200,bbox_inches='tight')
    plt.figure()
    samples=y.shape[1]
    ndim=y.shape[0]
    if (x is None):
        x = np.arange(1,ndim+1)
    if (usecolors):
        hsv = plt.get_cmap('hsv')
        colors = hsv(np.arange(0, 1, 1.0/samples))
    for i in range(samples):
        if (usecolors):
            color = colors[i]
        else:
            color='b'
        plt.plot(x,y[:,i],c=color,lw=lw)
        if (dots):
            plt.plot(x,y[:,i],'*',c=color,markersize=lw*3)
    plt.xlabel(xname)
    plt.ylabel(yname)
    if (xlim):
        plt.xlim(1,xlim)
    if (ylim):
        plt.ylim(-ylim,ylim)
    if (data is not None):
        plt.plot(data[0],data[1],'*',c='b',markersize=10)
    if (filename):
        plt.savefig(str(filename)+'.png', dpi=200,bbox_inches='tight')
    plt.figure()
    plt.errorbar(x,np.mean(y,axis=1),yerr=np.std(y,axis=1),ecolor='r')
    plt.plot(x,np.mean(y,axis=1),c='r')
    plt.xlabel(xname)
    plt.ylabel(yname)
    if (xlim):
        plt.xlim(1,xlim)
    if (ylim):
        plt.ylim(-ylim,ylim)
    if (filename):
        plt.savefig(str(filename)+'_erbr.png', dpi=200,bbox_inches='tight')


def inference(sample_x, sample_y, x, cov, ndim=40, length=6, sigma_v=0, l=2, sigma_f=1, sigma_n=0):
    nsample = sample_x.shape[0]
    k_data = np.zeros((ndim,nsample))
    cov_data = np.zeros((nsample,nsample))

    for i in range(ndim):
        for j in range(nsample):
            k_data[i,j]=k(sample_x[j],x[i],sigma_f,l)
    for i in range(nsample):
        for j in range(nsample):
            cov_data[i,j] = k(sample_x[i],sample_x[j],sigma_f,l) + sigma_v**2*direct_delta(x[i],x[j])
    mu = np.dot(np.dot(k_data,np.linalg.inv(cov_data)),sample_y)
    cov = cov - np.dot(np.dot(k_data,np.linalg.inv(cov_data)),k_data.T)
    return mu,cov

def k2(i,j,sigma_f,l):
    return np.power(sigma_f,2)*np.exp(-(np.power((i[0]-j[0]),2)/(2*np.power(l,2)) + np.power((i[1]-j[1]),2)/(2*np.power(l,2))))

def direct_delta2(i,j):
    if (i[0]==j[0] and i[1]==j[1]):
        return 1
    else:
        return 0

def sigma2(ndim=6,length=6,uniform=False,sigma_v=0,l=1,sigma_f=5):
    if (uniform):
        x = np.random.uniform(1,length,ndim)
        x.sort()
        y = np.random.uniform(1,length,ndim)
        y.sort()
    else:
        x = np.linspace(-length, length, ndim)
        y = np.linspace(-length, length, ndim)
    [X,Y] = np.meshgrid(x,y)
    X = X.flatten()
    Y = Y.flatten()
    cov = np.zeros((ndim*ndim,ndim*ndim))
    for i in range(ndim*ndim):
        for j in range(ndim*ndim):
                    cov[i,j]=k2((X[i],Y[i]),(X[j],Y[j]),sigma_f,l)+np.power(sigma_v,2)*direct_delta(i,j)
    return cov,X,Y

def inference2D(sample_x, sample_y, sample_z, x, y, cov, ndim=40, length=6, sigma_v=0, l=1, sigma_f=5, sigma_n=0):
    nsample = sample_x.shape[0]
    k_data = np.zeros((ndim**2,nsample))
    cov_data = np.zeros((nsample,nsample))
    for i in range(ndim**2):
        for j in range(nsample):
            k_data[i,j]=k2((sample_x[j],sample_y[j]),(x[i],y[i]),sigma_f,l)
    for i in range(nsample):
        for j in range(nsample):
            cov_data[i,j] = k2((sample_x[i],sample_y[i]),(sample_x[j],sample_y[j]),sigma_f,l) + sigma_v**2*direct_delta((sample_x[i],sample_x[j]),(sample_y[i],sample_y[j]))
    mu = np.dot(np.dot(k_data,np.linalg.inv(cov_data)),sample_z)
    cov = cov - np.dot(np.dot(k_data,np.linalg.inv(cov_data)),k_data.T)
    return mu,cov


def genGP(samples=10, ndim=40, length=6, sigma_v=0, l=2, sigma_f=1, uniform=True, dots=False, lw=3, xname="X", yname="Y", ylim=None, mean=0, usecolors=True, filename=None):
    samples=samples
    ndim=ndim
    length=length

    mu=np.ones(ndim)*mean
    cov,x=sigma(ndim,uniform=uniform,length=length,sigma_v=sigma_v,l=l,sigma_f=sigma_f)

    y=np.random.multivariate_normal(mu, cov, samples).T  
    plotGP(x=x,y=y,dots=dots,usecolors=usecolors,lw=lw,filename=filename,xname=xname,yname=yname,xlim=length,ylim=ylim)
    return x,y,cov

def closestReading(m,n,xs,ys,zs):
    xs = xs.flatten()
    ys = ys.flatten()
    zs = zs.flatten()
    aux = []
    for x, y in zip(xs,ys):
        aux.append(abs((m-x)**2 + (n-y)**2))
    return zs[aux.index(min(aux))]

def saveData(xs,ys,zs,filename):
    xs = xs.flatten()
    ys = ys.flatten()
    zs = zs.flatten()
    np.savetxt(filename, (xs,ys,zs), fmt='%1.4f')

def loadData(filename):
    return np.loadtxt(filename)

def debug():
    import pdb; pdb.set_trace()
