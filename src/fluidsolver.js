(function (window) {
    'use strict';

    // Export constructor
    window['FluidSolver'] = FluidSolver;

    /**
     * A Simple fluid solver implementation in javascript.
     *
     * Largely based on Jos Stam's paper "Real-Time Fluid Dynamics for Games".
     * @link http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
     *
     * Simulates the Navier–Stokes equations for incompressible fluids.
     * @link http://en.wikipedia.org/wiki/Navier-Stokes_equations
     *
     * Other implementations I've looked at while making this.
     * @link http://www.multires.caltech.edu/teaching/demos/java/stablefluids.htm
     * @link http://blog.inspirit.ru/fluidsolver-as3-port-of-msafluid/
     * @link http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html
     *
     * @author Topaz Bar <topaz1008@gmail.com>
     *
     * @param n {Number} Number of fluid cells for the simulation grid in each dimension (NxN)
     * @constructor
     */
    function FluidSolver(n) {
        var i;

        this.n = n;

        this.dt = 0.1; // The simulation time-step
        this.diffusion = 0.0002; // The amount of diffusion
        this.viscosity = 0; // The fluid's viscosity

        // Number of iterations to use in the Gauss-Seidel method in linearSolve()
        this.iterations = 10;

        this.doVorticityConfinement = true;
        this.doBuoyancy = true;

        // Two extra cells in each dimension for the boundaries
        this.numOfCells = (n + 2) * (n + 2) * (n + 2);

        this.tmp = null; // Scratch space for references swapping

        // This might benefit from using typed arrays like Float32Array in some configuration.
        // But I haven't seen any significant improvement on Chrome because V8 probably does it on its own.

        // Values for current simulation step
        this.u = new Array(this.numOfCells); // Velocity x
        this.v = new Array(this.numOfCells); // Velocity y
        this.w = new Array(this.numOfCells); // Velocity z
        this.d = new Array(this.numOfCells); // Density

        // Values from the last simulation step
        this.uOld = new Array(this.numOfCells);
        this.vOld = new Array(this.numOfCells);
        this.wOld = new Array(this.numOfCells);
        this.dOld = new Array(this.numOfCells);

        this.curlData = new Array(this.numOfCells); // The cell's curl

        // Initialize everything to zero
        for (i = 0; i < this.numOfCells; i++) {
            this.d[i] = this.u[i] = this.v[i] = this.w[i] = 0;
            this.dOld[i] = this.uOld[i] = this.vOld[i] = this.wOld[i] = 0;
            this.curlData[i] = 0;
        }
    }

    /** Boundaries enumeration. */
    FluidSolver.BOUNDARY_NONE = 0;
    FluidSolver.BOUNDARY_LEFT_RIGHT = 1;
    FluidSolver.BOUNDARY_TOP_BOTTOM = 2;
    FluidSolver.BOUNDARY_IN_OUT = 3;

    /**
     * A 'Private' stand alone function so closure compiler can inline this calculation once compiled.
     * DOES NOT make sure the indexes are integers.
     *
     * @param n {Number}
     * @param i {Number}
     * @param j {Number}
     * @param k {Number}
     * @returns {Number}
     * @private
     */
    function I(n, i, j, k) {
        return i + (n + 2) * j + (n + 2) * (n + 2) * k;
    }

    /**
     * Fluid cell indexing helper function.
     * (x | x) is a faster Math.floor(x)
     *
     * For public use.
     *
     * @return {number}
     * @public
     */
    FluidSolver.prototype.I = function (i, j, k) {
        return (i | i) + (this.n + 2) * (j | j) + (this.n + 2) * (this.n + 2) * (k | k);
    };

    /**
     * Density step.
     */
    FluidSolver.prototype.densityStep = function () {
        var i;

        this.addSource(this.d, this.dOld);

        this.swapD();
        this.diffuse(FluidSolver.BOUNDARY_NONE, this.d, this.dOld, this.diffusion);

        this.swapD();
        this.advect(FluidSolver.BOUNDARY_NONE, this.d, this.dOld, this.u, this.v, this.w);

        // Reset for next step
        for (i = 0; i < this.numOfCells; i++) {
            this.dOld[i] = 0;
        }
    };

    /**
     * Velocity step.
     */
    FluidSolver.prototype.velocityStep = function () {
        var i;

        this.addSource(this.u, this.uOld);
        this.addSource(this.v, this.vOld);
        this.addSource(this.w, this.wOld);

        if (this.doVorticityConfinement) {
            this.vorticityConfinement(this.uOld, this.vOld);
            this.addSource(this.u, this.uOld);
            this.addSource(this.v, this.vOld);
        }

        if (this.doBuoyancy) {
            this.buoyancy(this.vOld);
            this.addSource(this.v, this.vOld);
        }

        this.swapU();
        this.diffuse(FluidSolver.BOUNDARY_LEFT_RIGHT, this.u, this.uOld, this.viscosity);

        this.swapV();
        this.diffuse(FluidSolver.BOUNDARY_TOP_BOTTOM, this.v, this.vOld, this.viscosity);

        this.swapW();
        this.diffuse(FluidSolver.BOUNDARY_IN_OUT, this.w, this.wOld, this.viscosity);

        this.project(this.u, this.v, this.w, this.uOld, this.vOld);
        var a = this.v;
        this.swapU();
        this.swapV();
        this.swapW();

        this.advect(FluidSolver.BOUNDARY_LEFT_RIGHT, this.u, this.uOld, this.uOld, this.vOld, this.wOld);
        this.advect(FluidSolver.BOUNDARY_TOP_BOTTOM, this.v, this.vOld, this.uOld, this.vOld, this.wOld);
        this.advect(FluidSolver.BOUNDARY_IN_OUT, this.w, this.wOld, this.uOld, this.vOld, this.wOld);

        this.project(this.u, this.v, this.w, this.uOld, this.vOld);

        // Reset for next step
        for (i = 0; i < this.numOfCells; i++) {
            this.uOld[i] = this.vOld[i] = this.wOld[i] = 0;
        }
    };

    /**
     * Resets the density.
     */
    FluidSolver.prototype.resetDensity = function () {
        var i;
        for (i = 0; i < this.numOfCells; i++) {
            this.d[i] = 0;
        }
    };

    /**
     * Resets the velocity.
     */
    FluidSolver.prototype.resetVelocity = function () {
        var i;
        for (i = 0; i < this.numOfCells; i++) {
            // Set a small value so we can render the velocity field
            this.v[i] = this.u[i] = this.w[i] = 0.001;
        }
    };

    /**
     * Swap velocity x reference.
     * @private
     */
    FluidSolver.prototype.swapU = function () {
        this.tmp = this.u;
        this.u = this.uOld;
        this.uOld = this.tmp;
    };

    /**
     * Swap velocity y reference.
     * @private
     */
    FluidSolver.prototype.swapV = function () {
        this.tmp = this.v;
        this.v = this.vOld;
        this.vOld = this.tmp;
    };

    /**
     * Swap velocity z reference.
     * @private
     */
    FluidSolver.prototype.swapW = function () {
        this.tmp = this.w;
        this.w = this.wOld;
        this.wOld = this.tmp;
    };

    /**
     * Swap density reference.
     * @private
     */
    FluidSolver.prototype.swapD = function () {
        this.tmp = this.d;
        this.d = this.dOld;
        this.dOld = this.tmp;
    };

    /**
     * Integrate the density sources.
     *
     * @param x {Array<Number>}
     * @param s {Array<Number>}
     * @private
     */
    FluidSolver.prototype.addSource = function (x, s) {
        var i;
        for (i = 0; i < this.numOfCells; i++) {
            x[i] += s[i] * this.dt;
        }
    };

    /**
     * Calculate the curl at cell (i, j)
     * This represents the vortex strength at the cell.
     * Computed as: w = (del x U) where U is the velocity vector at (i, j).
     *
     * @param i Number
     * @param j {Number}
     * @return {Number}
     * @private
     */
    FluidSolver.prototype.curl = function (i, j) {
        var duDy = (this.u[I(this.n, i, j + 1)] - this.u[I(this.n, i, j - 1)]) * 0.5,
            dvDx = (this.v[I(this.n, i + 1, j)] - this.v[I(this.n, i - 1, j)]) * 0.5;

        return duDy - dvDx;
    };

    /**
     * Calculate the vorticity confinement force for each cell.
     * Fvc = (N x W) where W is the curl at (i, j) and N = del |W| / |del |W||.
     * N is the vector pointing to the vortex center, hence we
     * add force perpendicular to N.
     *
     * @param vcX {Array<Number>}
     * @param vcY {Array<Number>}
     * @private
     */
    FluidSolver.prototype.vorticityConfinement = function (vcX, vcY) {
        var i, j, dx, dy, norm, v;

        // Calculate magnitude of curl(i, j) for each cell
        for (i = 1; i <= this.n; i++) {
            for (j = 1; j <= this.n; j++) {
                this.curlData[I(this.n, i, j)] = Math.abs(this.curl(i, j));
            }
        }

        for (i = 1; i <= this.n; i++) {
            for (j = 1; j <= this.n; j++) {
                // Calculate the derivative of the magnitude (n = del |w|)
                dx = (this.curlData[I(this.n, i + 1, j)] - this.curlData[I(this.n, i - 1, j)]) * 0.5;
                dy = (this.curlData[I(this.n, i, j + 1)] - this.curlData[I(this.n, i, j - 1)]) * 0.5;

                norm = Math.sqrt((dx * dx) + (dy * dy));
                if (norm === 0) {
                    // Avoid divide by zero
                    norm = 1;
                }

                dx /= norm;
                dy /= norm;

                v = this.curl(i, j);

                // N x W
                vcX[I(this.n, i, j)] = dy * v * -1;
                vcY[I(this.n, i, j)] = dx * v;
            }
        }
    };

    /**
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
     */
    FluidSolver.prototype.buoyancy = function (buoy) {
        var i, j, k, length,
            tAmb = 0,
            a = 0.625,
            b = 0.25;

        // Sum all temperatures
//    for (i = 1; i <= this.n; i++) {
//        for (j = 1; j <= this.n; j++) {
//            tAmb += this.d[I(this.n, i, j)];
//        }
//    }

        // Sum all temperatures (faster)
        length = this.d.length;
        for (i = 0; i < length; i++) {
            tAmb += this.d[i];
        }

        // Calculate average temperature of the grid
        tAmb /= (this.n * this.n * this.n);

        // For each cell compute buoyancy force
        for (i = 1; i <= this.n; i++) {
            for (j = 1; j <= this.n; j++) {
                for (k = 1; k <= this.n; k++) {
                    buoy[I(this.n, i, j, k)] = a * this.d[I(this.n, i, j, k)] + -b * (this.d[I(this.n, i, j, k)] - tAmb);
                }
            }
        }
    };

    /**
     * Diffuse the density between neighbouring cells.
     *
     * @param b {Number}
     * @param x {Array<Number>}
     * @param x0 {Array<Number>}
     * @param diffusion {Number}
     * @private
     */
    FluidSolver.prototype.diffuse = function (b, x, x0, diffusion) {
        var a = this.dt * diffusion * this.n * this.n * this.n;

        this.linearSolve(b, x, x0, a, 1 + 6 * a);
    };

    /**
     * The advection step moves the density through the static velocity field.
     * Instead of moving the cells forward in time, we treat the cell's center as a particle
     * and then trace it back in time to look for the 'particles' which end up at the cell's center.
     *
     * @param b {Number}
     * @param d {Array<Number>}
     * @param d0 {Array<Number>}
     * @param u {Array<Number>}
     * @param v {Array<Number>}
     * @private
     */
    FluidSolver.prototype.advect = function (b, d, d0, u, v, w) {
        var i, j, k, i0, j0, k0, i1, j1, k1;
        var x, y, z, s0, t0, u0, s1, t1, u1, dt0;

        dt0 = this.dt * this.n;
        for (i = 1; i <= this.n; i++) {
            for (j = 1; j <= this.n; j++) {
                for (k = 1; k <= this.n; k++) {
                    x = i - dt0 * u[I(this.n, i, j, k)];
                    y = j - dt0 * v[I(this.n, i, j, k)];
                    z = k - dt0 * w[I(this.n, i, j, k)];

                    if (x < 0.5) x = 0.5;
                    if (x > this.n + 0.5) x = this.n + 0.5;

                    i0 = (x | x);
                    i1 = i0 + 1;

                    if (y < 0.5) y = 0.5;
                    if (y > this.n + 0.5) y = this.n + 0.5;

                    j0 = (y | y);
                    j1 = j0 + 1;

                    if (z < 0.5) z = 0.5;
                    if (z > this.n + 0.5) z = this.n + 0.5;

                    k0 = (z | z);
                    k1 = k0 + 1;

                    s1 = x - i0;
                    s0 = 1 - s1;
                    t1 = y - j0;
                    t0 = 1 - t1;
                    u1 = z - k0;
                    u0 = 1 - u1;

                    d[I(this.n, i, j, k)] = s0 * (t0 * u0 * d0[I(this.n, i0, j0, k0)] + t1 * u0 * d0[I(this.n, i0, j1, k0)]) +
                                            s1 * (t0 * u0 * d0[I(this.n, i1, j0, k0)] + t1 * u0 * d0[I(this.n, i1, j1, k0)]) +
                                            s0 * (t0 * u1 * d0[I(this.n, i0, j0, k1)] + t1 * u1 * d0[I(this.n, i0, j1, k1)]) +
                                            s1 * (t0 * u1 * d0[I(this.n, i1, j0, k1)] + t1 * u1 * d0[I(this.n, i1, j1, k1)]);

                }
            }
        }

        this.setBoundary(b, d);
    };

    /**
     * Forces the velocity field to be mass conserving.
     * This step is what actually produces the nice looking swirly vortices.
     *
     * It uses a result called Hodge Decomposition which says that every velocity field is the sum
     * of a mass conserving field, and a gradient field. So we calculate the gradient field, and subtract
     * it from the velocity field to get a mass conserving one.
     * It solves a linear system of equations called Poisson Equation.
     *
     * @param u {Array<Number>}
     * @param v {Array<Number>}
     * @param w {Array<Number>}
     * @param p {Array<Number>}
     * @param div {Array<Number>}
     * @private
     */
    FluidSolver.prototype.project = function (u, v, w, p, div) {
        var i, j, k;

        // Calculate the gradient field
        var h = 1.0 / this.n;
        for (i = 1; i <= this.n; i++) {
            for (j = 1; j <= this.n; j++) {
                for (k = 1; k <= this.n; k++) {
                    div[I(this.n, i, j, k)] = -1.0/3.0 * h * (u[I(this.n, i + 1, j, k)] - u[I(this.n, i - 1, j, k)] + v[I(this.n, i, j + 1, k)] - v[I(this.n, i, j - 1, k)] + w[I(this.n, i, j, k + 1)] - w[I(this.n, i, j, k - 1)]);
                    p[I(this.n, i, j, k)] = 0;
                }
            }
        }

        this.setBoundary(FluidSolver.BOUNDARY_NONE, div);
        this.setBoundary(FluidSolver.BOUNDARY_NONE, p);

        // Solve the Poisson equations
        this.linearSolve(FluidSolver.BOUNDARY_NONE, p, div, 1, 6);

        // Subtract the gradient field from the velocity field to get a mass conserving velocity field.
        for (i = 1; i <= this.n; i++) {
            for (j = 1; j <= this.n; j++) {
                for (k = 1; k <= this.n; k++) {
                    u[I(this.n, i, j, k)] -= 0.5 * (p[I(this.n, i + 1, j, k)] - p[I(this.n, i - 1, j, k)]) / h;
                    v[I(this.n, i, j, k)] -= 0.5 * (p[I(this.n, i, j + 1, k)] - p[I(this.n, i, j - 1, k)]) / h;
                    w[I(this.n, i, j, k)] -= 0.5 * (p[I(this.n, i, j, k + 1)] - p[I(this.n, i, j, k - 1)]) / h;
                }
            }
        }

        this.setBoundary(FluidSolver.BOUNDARY_LEFT_RIGHT, u);
        this.setBoundary(FluidSolver.BOUNDARY_TOP_BOTTOM, v);
        this.setBoundary(FluidSolver.BOUNDARY_IN_OUT, w);
    };

    /**
     * Solve a linear system of equations using Gauss-Seidel method.
     *
     * @param b {Number}
     * @param x {Array<Number>}
     * @param x0 {Array<Number>}
     * @param a {Number}
     * @param c {Number}
     * @private
     */
    FluidSolver.prototype.linearSolve = function (b, x, x0, a, c) {
        var i, j, k, l, invC = 1.0 / c;

        for (l = 0; l < this.iterations; l++) {
            for (i = 1; i <= this.n; i++) {
                for (j = 1; j <= this.n; j++) {
                    for (k = 1; k <= this.n; k++) {
                        x[I(this.n, i, j, k)] = (x0[I(this.n, i, j, k)] + a * (x[I(this.n, i - 1, j, k)] + x[I(this.n, i + 1, j, k)] +
                                                                               x[I(this.n, i, j - 1, k)] + x[I(this.n, i, j + 1, k)] +
                                                                               x[I(this.n, i, j, k - 1)] + x[I(this.n, i, j, k + 1)])) * invC;
                    }
                }
            }

            this.setBoundary(b, x);
        }
    };

    /**
     * Set boundary conditions.
     *
     * @param b {Number}
     * @param x {Array<Number>}
     * @private
     */
    FluidSolver.prototype.setBoundary = function (b, x) {
        var i, j;

        for (i = 1; i <= this.n; i++) {
            for (j = 1; j <= this.n; j++) {

                x[I(this.n, 0, i, j)] = (b === FluidSolver.BOUNDARY_LEFT_RIGHT) ? -x[I(this.n, 1, i, j)] : x[I(this.n, 1, i, j)];

                x[I(this.n, this.n + 1, i, j)] = (b === FluidSolver.BOUNDARY_LEFT_RIGHT) ? -x[I(this.n, this.n, i, j)] : x[I(this.n, this.n, i, j)];

                x[I(this.n, i, 0, j)] = (b === FluidSolver.BOUNDARY_TOP_BOTTOM) ? -x[I(this.n, i, 1, j)] : x[I(this.n, i, 1, j)];

                x[I(this.n, i, this.n + 1, j)] = (b === FluidSolver.BOUNDARY_TOP_BOTTOM) ? -x[I(this.n, i, this.n, j)] : x[I(this.n, i, this.n, j)];

                x[I(this.n, i, j, 0)] = (b === FluidSolver.BOUNDARY_IN_OUT) ? -x[I(this.n, i, j, 1)] : x[I(this.n, i, j, 1)];

                x[I(this.n, i, j, this.n + 1)] = (b === FluidSolver.BOUNDARY_IN_OUT) ? -x[I(this.n, i, j, this.n)] : x[I(this.n, i, j, this.n)];
            }
        }

        x[I(this.n, 0, 0, 0)] = 1.0/3.0 * (x[I(this.n, 1, 0, 0)] + x[I(this.n, 0, 1, 0)] + x[I(this.n, 0, 0, 1)]);
        x[I(this.n, 0, this.n + 1, 0)] = 1.0/3.0 * (x[I(this.n, 1, this.n + 1, 0)] + x[I(this.n, 0, this.n, 0)] + x[I(this.n, 0, this.n + 1, 1)]);
        x[I(this.n, this.n + 1, 0, 0)] = 1.0/3.0 * (x[I(this.n, this.n, 0, 0)] + x[I(this.n, this.n + 1, 1, 0)] + x[I(this.n, this.n + 1, 0, 1)]);
        x[I(this.n, this.n + 1, this.n + 1, 0)] = 1.0/3.0 * (x[I(this.n, this.n, this.n + 1, 0)] + x[I(this.n, this.n + 1, this.n, 0)] + x[I(this.n, this.n + 1, this.n + 1, 1)]);
        x[I(this.n, 0, 0, this.n + 1)] = 1.0/3.0 * (x[I(this.n, 1, 0, this.n + 1)] + x[I(this.n, 0, 1, this.n + 1)] + x[I(this.n, 0, 0, this.n)]);
        x[I(this.n, 0, this.n + 1, this.n + 1)] = 1.0/3.0 * (x[I(this.n, 1, this.n + 1, this.n + 1)] + x[I(this.n, 0, this.n, this.n + 1)] + x[I(this.n, 0, this.n + 1, this.n)]);
        x[I(this.n, this.n + 1, 0, this.n + 1)] = 1.0/3.0 * (x[I(this.n, this.n, 0, this.n + 1)] + x[I(this.n, this.n + 1, 1, this.n + 1)] + x[I(this.n, this.n + 1, 0, this.n)]);
        x[I(this.n, this.n + 1, this.n + 1, this.n + 1)] = 1.0/3.0 * (x[I(this.n, this.n, this.n + 1, this.n + 1)] + x[I(this.n, this.n + 1, this.n, this.n + 1)] + x[I(this.n, this.n + 1, this.n + 1, this.n)]);




    };

})(window);
