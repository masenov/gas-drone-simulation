import numpy as np
import time
import timeit
import matplotlib.pyplot as plt
import h5py
import random

class FluidSimulation:
    def __init__(self, n, windDirection=180.0, windLocations=32, windSpeed=0.02, windNoise=180, windNoiseTimestep=300, windLocalNoise=0.001):

        self.n = n

        self.dt = 0.1 # The simulation time-step
        self.diffusion = 0.0001 # The amount of diffusion
        self.viscosity = 0 # The fluid's viscosity

        # Number of iterations to use in the Gauss-Seidel method in linearSolve()
        self.iterations = 10

        self.doVorticityConfinement = False
        self.doBuoyancy = False

        # Two extra cells in each dimension for the boundaries
        self.numOfCells = (n + 2) * (n + 2)

        self.tmp = None # Scratch space for references swapping

        self.windSpeed = windSpeed
        self.windDirection = windDirection
        self.windLocations = windLocations
        # Change wind direction after N-th timesteps
        self.windNoise = windNoise
        self.windNoiseTimestep = windNoiseTimestep
        self.windNoiseCurrentStep = self.windNoiseTimestep
        self.windNoiseCurrent = 0
        # Local noise at each location
        self.windLocalNoise = windLocalNoise

        self.gasRelease = 110
        self.gasLocationX = n/2
        self.gasLocationY = n/2
        # Self might benefit from using typed arrays like Float32Array in some configuration.
        # But I haven't seen any significant improvement on Chrome because V8 probably does it on its own.

        # Values for current simulation step
        self.u = np.ones(self.numOfCells) * 0.001 # Velocity x
        self.v = np.ones(self.numOfCells) * 0.001 # Velocity y
        self.d = np.zeros(self.numOfCells) # Density

        # Values from the last simulation step
        self.uOld = np.zeros(self.numOfCells)
        self.vOld = np.zeros(self.numOfCells)
        self.dOld = np.zeros(self.numOfCells)

        self.curlData = np.zeros(self.numOfCells) # The cell's curl

        # Initialize everything to zero
        for i in range(self.numOfCells):
            self.d[i] = self.u[i] = self.v[i] = 0
            self.dOld[i] = self.uOld[i] = self.vOld[i] = 0
            self.curlData[i] = 0

        # Boundaries enumeration
        self.BOUNDARY_NONE = 0
        self.BOUNDARY_LEFT_RIGHT = 1
        self.BOUNDARY_TOP_BOTTOM = 2
        self.NUM_OF_CELLS = n # Number of cells (not including the boundary)
        self.VIEW_SIZE = 640    # View size (square)
        self.FPS = 60           # Frames per second


        self.CELL_SIZE = self.VIEW_SIZE / self.NUM_OF_CELLS  # Size of each cell in pixels
        self.CELL_SIZE_CEIL = np.ceil(self.CELL_SIZE) # Size of each cell in pixels (ceiling)

        self.i = np.arange(1,self.n+1)
        self.j = np.arange(1,self.n+1)
        ii, jj = np.meshgrid(self.i,self.j)
        self.ii = ii.ravel()
        self.jj = jj.ravel()
        self.Iij = self.I(self.ii, self.jj)
        print (self.ii.shape, self.Iij.shape)
        self.Iij0 = self.I(self.ii + 1, self.jj)
        self.Iij1 = self.I(self.ii - 1, self.jj)
        self.Iij2 = self.I(self.ii, self.jj + 1)
        self.Iij3 = self.I(self.ii, self.jj - 1)
        self.index = 0

    def randomWind(self):
        return (random.random()*2-1)*self.windLocalNoise

    def I(self, i, j):
        return i + (self.n + 2) * j

    """
     * Density step.
     """
    def densityStep(self):

        self.addSource(self.d, self.dOld)

        self.swapD()
        self.diffuse(self.BOUNDARY_NONE, self.d, self.dOld, self.diffusion)

        self.swapD()
        self.advect(self.BOUNDARY_NONE, self.d, self.dOld, self.u, self.v)

        # Reset for next step
        self.dOld.fill(0)

    """
     * Velocity step.
     """
    def velocityStep(self):
        self.addSource(self.u, self.uOld)
        self.addSource(self.v, self.vOld)

        if (self.doVorticityConfinement):
            self.vorticityConfinement(self.uOld, self.vOld)
            self.addSource(self.u, self.uOld)
            self.addSource(self.v, self.vOld)

        if (self.doBuoyancy):
            self.buoyancy(self.vOld)
            self.addSource(self.v, self.vOld)
        self.swapU()
        self.diffuse(self.BOUNDARY_LEFT_RIGHT, self.u, self.uOld, self.viscosity)

        self.swapV()
        self.diffuse(self.BOUNDARY_TOP_BOTTOM, self.v, self.vOld, self.viscosity)

        self.project(self.u, self.v, self.uOld, self.vOld)
        self.swapU()
        self.swapV()

        self.advect(self.BOUNDARY_LEFT_RIGHT, self.u, self.uOld, self.uOld, self.vOld)
        self.advect(self.BOUNDARY_TOP_BOTTOM, self.v, self.vOld, self.uOld, self.vOld)

        self.project(self.u, self.v, self.uOld, self.vOld)

        # Reset for next step
        self.uOld.fill(0)
        self.vOld.fill(0)

    """
     * Resets the density.
     """
    def resetDensity(self):
        self.d.fill(0)

    """
     * Resets the velocity.
     """
    def resetVelocity(self):
        self.v.fill(0.001)
        self.u.fill(0.001)

    """
     * Swap velocity x reference.
     """
    def swapU(self):
        self.u, self.uOld = self.uOld, self.u

    """
     * Swap velocity y reference.
     """
    def swapV(self):
        self.v, self.vOld = self.vOld, self.v

    """
     * Swap density reference.
     """
    def swapD(self):
        self.d, self.dOld = self.dOld, self.d

    """
     * Integrate the density sources.
     """
    def addSource(self, x, s):
        x += s * self.dt

    """
     * Calculate the curl at cell (i, j)
     * This represents the vortex strength at the cell.
     * Computed as: w = (del x U) where U is the velocity vector at (i, j).
     """
    def curl(self, i, j):
        duDy = (self.u[self.I(i, j + 1)] - self.u[self.I(i, j - 1)]) * 0.5
        dvDx = (self.v[self.I(i + 1, j)] - self.v[self.I(i - 1, j)]) * 0.5

        return duDy - dvDx

    """
     * Calculate the vorticity confinement force for each cell.
     * Fvc = (N x W) where W is the curl at (i, j) and N = del |W| / |del |W||.
     * N is the vector pointing to the vortex center, hence we
     * add force perpendicular to N.
     """
    def vorticityConfinement(self, vcX, vcY):

        # Calculate magnitude of curl(i, j) for each cell
        for i in range(self.n):
            for j in range(self.n):
                self.curlData[self.I(i, j)] = np.abs(self.curl(i, j))

        for i in range(self.n):
            for j in range(self.n):
                # Calculate the derivative of the magnitude (n = del |w|)
                dx = (self.curlData[self.I(i + 1, j)] - self.curlData[self.I(i - 1, j)]) * 0.5
                dy = (self.curlData[self.I(i, j + 1)] - self.curlData[self.I(i, j - 1)]) * 0.5

                norm = np.sqrt((dx * dx) + (dy * dy))
                if (norm == 0):
                    # Avoid divide by zero
                    norm = 1

                dx /= norm
                dy /= norm

                v = self.curl(i, j)

                # N x W
                vcX[self.I(i, j)] = dy * v * -1
                vcY[self.I(i, j)] = dx * v

    """
     * Calculate the buoyancy force for the grid.
     * Fbuoy = -a * d * Y + b * (T - Tamb) * Y where Y = (0,1)
     * The constants a and b are positive with physically meaningful quantities.
     * T is the temperature at the current cell, Tamb is the average temperature of the fluid grid
     *
     * In this simplified implementation we say that the temperature is synonymous with density
     * and because there are no other heat sources we can just use the density field instead of adding a new
     * temperature field.
     *
     * @param buoy {Array<Number>}
     * @private
     """
    def buoyancy(self, buoy):
        tAmb = 0
        a = 0.000625
        b = 0.025

        # Sum all temperatures
#    for i in range(self.n):
#        for j in range(self.n):
#            tAmb += self.d[I(self.n, i, j)]
#        }
#    }

        # Sum all temperatures (faster)
        length = self.d.length
        for i in range(length):
            tAmb += self.d[i]

        # Calculate average temperature of the grid
        tAmb /= (self.n * self.n)

        # For each cell compute buoyancy force
        for i in range(self.n):
            for j in range(self.n):
                buoy[I(self.n, i, j)] = a * self.d[I(self.n, i, j)] + -b * (self.d[I(self.n, i, j)] - tAmb)

    """
     * Diffuse the density between neighbouring cells.
     """
    def diffuse(self, b, x, x0, diffusion):
        a = self.dt * diffusion * self.n * self.n

        self.linearSolve(b, x, x0, a, 1 + 4 * a)

    """
     * The advection step moves the density through the static velocity field.
     * Instead of moving the cells forward in time, we treat the cell's center as a particle
     * and then trace it back in time to look for the 'particles' which end up at the cell's center.
     """
    def advect(self, b, d, d0, u, v):


        dt0 = self.dt * self.n
        """
        i = np.arange(1,self.n+1)
        j = np.arange(1,self.n+1)
        ii, jj = np.meshgrid(i,j)
        ii = ii.flatten()
        jj = jj.flatten()
        Iij = self.I(ii, jj)
        """
        x = self.ii - dt0 * u[self.Iij]
        y = self.jj - dt0 * v[self.Iij]
        x[x < 0.5] = 0.5
        x[x > self.n + 0.5] = self.n + 0.5
        i0 = x.astype(int)
        i1 = i0 + 1
        y[y < 0.5] = 0.5
        y[y > self.n + 0.5] = self.n + 0.5
        j0 = y.astype(int)
        j1 = j0 + 1
        s1 = x - i0
        s0 = 1 - s1
        t1 = y - j0
        t0 = 1 - t1
        Iij0 = self.I(i0, j0)
        Iij1 = self.I(i0, j1)
        Iij2 = self.I(i1, j0)
        Iij3 = self.I(i1, j1)
        d[self.Iij] = s0 * (t0 * d0[Iij0] + t1 * d0[Iij1]) + s1 * (t0 * d0[Iij2] + t1 * d0[Iij3])
        #print (ii.shape)
        """
        for i in range(self.n):
            for j in range(self.n):
                x = i - dt0 * u[self.I(i, j)]
                y = j - dt0 * v[self.I(i, j)]

                if (x < 0.5):
                    x = 0.5
                if (x > self.n + 0.5):
                    x = self.n + 0.5

                i0 = int(x)
                i1 = i0 + 1

                if (y < 0.5):
                    y = 0.5
                if (y > self.n + 0.5):
                    y = self.n + 0.5

                j0 = int(y)
                j1 = j0 + 1
                s1 = x - i0
                s0 = 1 - s1
                t1 = y - j0
                t0 = 1 - t1
                d[self.I(i, j)] = s0 * (t0 * d0[self.I(i0, j0)] + t1 * d0[self.I(i0, j1)]) + s1 * (t0 * d0[self.I(i1, j0)] + t1 * d0[self.I(i1, j1)])
        """
        self.setBoundary(b, d)

    """
     * Forces the velocity field to be mass conserving.
     * This step is what actually produces the nice looking swirly vortices.
     *
     * It uses a result called Hodge Decomposition which says that every velocity field is the sum
     * of a mass conserving field, and a gradient field. So we calculate the gradient field, and subtract
     * it from the velocity field to get a mass conserving one.
     * It solves a linear system of equations called Poisson Equation.
     *
     * @param u:Array<Number>}
     * @param v:Array<Number>}
     * @param p:Array<Number>}
     * @param div:Array<Number>}
     * @private
     """
    def project(self, u, v, p, div):

        # Calculate the gradient field
        h = 1.0 / self.n
        """
        i = np.arange(1,self.n+1)
        j = np.arange(1,self.n+1)
        ii, jj = np.meshgrid(i,j)
        ii = ii.flatten()
        jj = jj.flatten()
        Iij = self.I(ii, jj)
        Iij0 = self.I(ii + 1, jj)
        Iij1 = self.I(ii - 1, jj)
        Iij2 = self.I(ii, jj + 1)
        Iij3 = self.I(ii, jj - 1)
        """
        div[self.Iij] = -0.5 * h * (u[self.Iij0] - u[self.Iij1] + v[self.Iij2] - v[self.Iij3])
        p.fill(0.0)

        """
        for i in range(self.n):
            for j in range(self.n):
                div[self.I(i, j)] = -0.5 * h * (u[self.I(i + 1, j)] - u[self.I(i - 1, j)] +
                    v[self.I(i, j + 1)] - v[self.I(i, j - 1)])

                p[self.I(i, j)] = 0
        """
        self.setBoundary(self.BOUNDARY_NONE, div)
        self.setBoundary(self.BOUNDARY_NONE, p)

        # Solve the Poisson equations
        self.linearSolve(self.BOUNDARY_NONE, p, div, 1, 4)

        # Subtract the gradient field from the velocity field to get a mass conserving velocity field.


        u[self.Iij] -= 0.5 * (p[self.Iij0] - p[self.Iij1]) / h
        v[self.Iij] -= 0.5 * (p[self.Iij2] - p[self.Iij3]) / h
        """
        for i in range(self.n):
            for j in range(self.n):
                u[self.I(i, j)] -= 0.5 * (p[self.I(i + 1, j)] - p[self.I(i - 1, j)]) / h
                v[self.I(i, j)] -= 0.5 * (p[self.I(i, j + 1)] - p[self.I(i, j - 1)]) / h
        """
        self.setBoundary(self.BOUNDARY_LEFT_RIGHT, u)
        self.setBoundary(self.BOUNDARY_TOP_BOTTOM, v)

    """
     * Solve a linear system of equations using Gauss-Seidel method.
     *
     * @param b:Number}
     * @param x:Array<Number>}
     * @param x0:Array<Number>}
     * @param a:Number}
     * @param c:Number}
     * @private
     """
    def linearSolve(self, b, x, x0, a, c):
        invC = 1.0 / c
        for k in range(self.iterations):
            """
            i = np.arange(1,self.n+1)
            j = np.arange(1,self.n+1)
            ii, jj = np.meshgrid(i,j)
            ii = ii.flatten()
            jj = jj.flatten()
            Iij = self.I(ii, jj)
            Iij0 = self.I(ii - 1, jj)
            Iij1 = self.I(ii + 1, jj)
            Iij2 = self.I(ii, jj - 1)
            Iij3 = self.I(ii, jj + 1)
            """
            #print (ii.shape,jj.shape,Iij.shape)
            x[self.Iij] = (x0[self.Iij] + a*(x[self.Iij1] + x[self.Iij0] + x[self.Iij3] + x[self.Iij2])) * invC
            """
            for i in range(self.n):
                for j in range(self.n):
                    x[self.I(i, j)] = (x0[self.I(i, j)] + a * (x[self.I(i - 1, j)] + x[self.I(i + 1, j)] +
                        x[self.I(i, j - 1)] + x[self.I(i, j + 1)])) * invC
            """
            self.setBoundary(b, x)

    """
     * Set boundary conditions.
     *
     * @param b:Number}
     * @param x:Array<Number>}
     * @private
     """
    def setBoundary(self, b, x):


        #i = np.arange(1,self.n+1)

        if (b == self.BOUNDARY_LEFT_RIGHT):
            x[self.I(0, self.i)] = -x[self.I(1, self.i)]
        else:
            x[self.I(0, self.i)].fill(0.0)

        if (b == self.BOUNDARY_LEFT_RIGHT):
            x[self.I(self.n + 1, self.i)] = -x[self.I(self.n, self.i)]
        else:
            x[self.I(self.n + 1, self.i)].fill(0.0)

        if (b == self.BOUNDARY_TOP_BOTTOM):
            x[self.I(self.i, 0)] = -x[self.I(self.i, 1)]
        else:
            x[self.I(self.i, 0)].fill(0.0)

        if (b == self.BOUNDARY_TOP_BOTTOM):
            x[self.I(self.i, self.n + 1)] = -x[self.I(self.i, self.n)]
        else:
            x[self.I(self.i, self.n + 1)].fill(0.0)

        """
        for i in range(self.n):
            x[self.I(0, i)] = (b == self.BOUNDARY_LEFT_RIGHT) and -x[self.I(1, i)] or 0

            x[self.I(self.n + 1, i)] = (b == self.BOUNDARY_LEFT_RIGHT) and -x[self.I(self.n, i)] or 0

            x[self.I(i, 0)] = (b == self.BOUNDARY_TOP_BOTTOM) and -x[self.I(i, 1)] or 0

            x[self.I(i, self.n + 1)] = (b == self.BOUNDARY_TOP_BOTTOM) and -x[self.I(i, self.n)] or 0
        """
        x[self.I(0, 0)] = 0.5 * (x[self.I(1, 0)] + x[self.I(0, 1)])
        x[self.I(0, self.n + 1)] = 0.5 * (x[self.I(1, self.n + 1)] + x[self.I(0, self.n)])
        x[self.I(self.n + 1, 0)] = 0.5 * (x[self.I(self.n, 0)] + x[self.I(self.n + 1, 1)])
        x[self.I(self.n + 1, self.n + 1)] = 0.5 * (x[self.I(self.n, self.n + 1)] + x[self.I(self.n + 1, self.n)])

    def toRadians (self, angle):
        return angle * (np.pi/ 180)

    def update(self, f):
        invMaxColor = 1.0 / 255

        self.dOld[self.I(self.gasLocationX, self.gasLocationY)] = self.gasRelease
        print (self.d[self.I(10,10)])
        start = time.time()

        # Step the fluid simulation
        #print ("a")
        self.velocityStep()
        end = time.time()
       #print (end - start)
        #print ("b")
        self.densityStep()
        end2 = time.time()
       #print (end2 - end)
        #print ("c")



        # Artificial wind field
        if self.windNoiseCurrentStep == self.windNoiseTimestep:
            rand = (np.random.rand()*2 - 1)*self.windNoise
            self.windNoiseCurrent = rand
            self.windNoiseCurrentStep = 0
        else:
            rand = self.windNoiseCurrent
            self.windNoiseCurrentStep += 1
        du = np.sin(self.toRadians(self.windDirection + rand - 180))*self.windSpeed
        dv = np.cos(self.toRadians(self.windDirection + rand - 180))*self.windSpeed
        # Add the mouse velocity to cells above, below, to the left, and to the right as well.
        nps = self.windLocations
        acc = self.NUM_OF_CELLS/nps
        for i in range(nps):
            for j in range(nps):
                ii = int(i*acc)
                jj = int(j*acc)
                self.uOld[self.I(ii, jj)] = du + self.randomWind()
                self.vOld[self.I(ii, jj)] = dv + self.randomWind()
        #self.update()
        # lastTime is now
        # End update()
        #np.save('data/test' + str(self.index) + '.npy',self.d)
        dset = f.create_dataset('test' + str(self.index), (self.n + 2,self.n + 2), dtype='f', compression="gzip")
        dset[...] = self.d.reshape(self.n+2,self.n+2)
        dset2 = f.create_dataset('wind' + str(self.index), (1,2), dtype='f', compression="gzip")
        dset2[...] = np.array([[self.windDirection + rand, self.windSpeed]])
        self.index += 1


    def getRandom(self, min, max):
        return min + np.random() * ((max + 1) - min)


