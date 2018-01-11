/**
 * Demo usage of the FluidSolver class.
 *
 * @author Topaz Bar <topaz1008@gmail.com>
 */
(function (window, document) {
    'use strict';

    // setInterval still seems to be faster than this most of the time.
    window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame ||
        window.webkitRequestAnimationFrame || window.msRequestAnimationFrame;

    updateSimulation();
//    setInterval(update, 1000 / FPS);
    var NUM_OF_CELLS = 256, // Number of cells (not including the boundary)
        VIEW_SIZE = 640,    // View size (square)
        FPS = 60;           // Frames per second
    var timestep = 0;


    var CELL_SIZE = VIEW_SIZE / NUM_OF_CELLS,  // Size of each cell in pixels
        CELL_SIZE_CEIL = Math.ceil(CELL_SIZE); // Size of each cell in pixels (ceiling)

    // Globals
    var canvas = document.getElementById('main'),
        context = canvas.getContext('2d');

    // We draw the density on a bitmap for performance reasons
    var fdBuffer = context.createImageData(VIEW_SIZE, VIEW_SIZE);

    // Set render states
    canvas.width = canvas.height = VIEW_SIZE;       // View size
    context.lineWidth = 1;                          // Velocity field line width
    context.strokeStyle = 'rgb(192, 0, 0)';         // Velocity field color
    //context.globalCompositeOperation = 'screen';  // Blend mode

    // Disable smoothing when using floating point pixel values
    context.imageSmoothingEnabled = false;
    context.webkitImageSmoothingEnabled = false;
    context.mozImageSmoothingEnabled = false;


    // Demo app variables
    var lastTime = Date.now(),
        isMouseDown = false,
        oldMouseX = 0,
        oldMouseY = 0,
        particles = [];



    var options = {
        drawVelocityField: false,
        drawDensityField: true,
        drawParticles: true,
        grayscale: false,
        resetParticles: function () { particles.length = 0; }
    };


    /**
     * Main mouse move listener
     *
     * @param event {MouseEvent|Object}
     */
    /**
     * Draw the simulation grid.
     */
    function drawGrid() {
        var i;

        context.lineWidth = 1;
        context.strokeStyle = 'rgb(255, 255, 255)';
        context.beginPath();

        // Vertical
        for (i = 0; i <= VIEW_SIZE; i += CELL_SIZE_CEIL) {
            context.moveTo(i, 0);
            context.lineTo(i, VIEW_SIZE);
        }

        // Horizontal
        for (i = 0; i <= VIEW_SIZE; i += CELL_SIZE_CEIL) {
            context.moveTo(0, i);
            context.lineTo(VIEW_SIZE, i);
        }

        context.stroke();
    }

    /**
     * Clears all the pixels on the image data.
     *
     * @param image {ImageData}
     */
    function clearImageData(image) {
        var i, length = image.data.length;
        for (i = 0; i < length; i++) {
            image.data[i] = 0;
        }
    }

    function updateDrawing() {
        // Clear the canvas
        context.clearRect(0, 0, VIEW_SIZE, VIEW_SIZE);

        // Draw the last frame's buffer and clear for drawing the current.
        context.putImageData(fdBuffer, 0, 0);
        clearImageData(fdBuffer);

        //drawGrid();

        if (options.drawVelocityField) {
            // Call once per frame
            context.lineWidth = 1;
            context.strokeStyle = 'rgb(192, 0, 0)';
            context.beginPath();
        }


    }

    function updateSimulation() {
        var socket = io();
        socket.emit('update simulation', timestep);
        //alert('Send request');
        socket.on('receive simulation data', function(newfs){
            //alert("Received data");
            updateRender(newfs);
            socket.emit('update simulation', timestep);
        });
    }

    function updateRender(fs) {
        var i, j, k, l, m, dx, dy, color, cellIndex, deltaTime,
            r, g, b, u, v, density, pxIdx, pxX, pxY,
            invMaxColor = 1.0 / 255;
        // Render fluid
        for (i = 1; i <= NUM_OF_CELLS; i++) {
            // The x position of current cell
            dx = (i - 0.5) * CELL_SIZE;

            for (j = 1; j <= NUM_OF_CELLS; j++) {
                // The y position of current cell
                dy = (j - 0.5) * CELL_SIZE;

                cellIndex = i + (NUM_OF_CELLS + 2) * j;

                // Draw density
                density = fs.d[cellIndex];
                if (options.drawDensityField && density > 0) {
                    color = density * 255;

                    // fdBuffer.data is actually a Uint8ClampedArray so there is no need to manually clamp color values
                    //if (color < 0) color = 0;
                    //if (color > 255) color = 255;

                    r = color;
                    g = 20;
                    b = 20;

                    // Draw the cell on an image for performance reasons
                    for (l = 0; l < CELL_SIZE_CEIL; l++) {
                        for (m = 0; m < CELL_SIZE_CEIL; m++) {
                            pxX = (i - 1) * CELL_SIZE + l;
                            pxY = (j - 1) * CELL_SIZE + m;
                            pxIdx = ((pxX | pxX) + (pxY | pxY) * VIEW_SIZE) * 4;

                            fdBuffer.data[pxIdx    ] = r;
                            fdBuffer.data[pxIdx + 1] = g;
                            fdBuffer.data[pxIdx + 2] = b;
                            fdBuffer.data[pxIdx + 3] = 255;
                        }
                    }
                }

                // Draw velocity field ?
                if (options.drawVelocityField && (i % 2) === 0 && (j % 2) === 0) {
                    u = fs.u[cellIndex] * 50;
                    v = fs.v[cellIndex] * 50;

                    context.moveTo(dx, dy);
                    context.lineTo(dx + u, dy + v);
                }

            } // End for all cells in the y direction

        } // End for all cells in the x direction

        if (options.drawVelocityField) {
            // Call once per frame
            context.stroke();
        }

        // Update and render particles
        var x0, y0, p, alpha, lastAlpha = 0, particlesLength = particles.length;
        if (options.drawParticles) {
            context.lineWidth = 2;
            context.strokeStyle = 'rgb(255, 255, 255)';
            context.beginPath();

            for (k = 0; k < particlesLength; k++) {
                p = particles[k];

                p.age += deltaTime;
                alpha = (1 - (p.age / Particle.TIME_TO_LIVE));
                if ((alpha < 0.01) ||
                    (p.age >= Particle.TIME_TO_LIVE) ||
                    (p.x <= 0 || p.x >= VIEW_SIZE || p.y <= 0 || p.y >= VIEW_SIZE)) {
                    p.dead = true;

                } else {
                    x0 = (p.x / VIEW_SIZE) * NUM_OF_CELLS + 2;
                    y0 = (p.y / VIEW_SIZE) * NUM_OF_CELLS + 2;

                    cellIndex = fs.I(x0, y0);

                    p.vx = fs.u[cellIndex] * 50;
                    p.vy = fs.v[cellIndex] * 50;

                    p.x += p.vx;
                    p.y += p.vy;

                    if (Math.abs(alpha - lastAlpha) > 0.001) {
                        // Only change stroke style if the alpha changed to save on render state changes.
                        context.strokeStyle = 'rgba(255, 255, 255, ' + alpha + ')';
                        lastAlpha = alpha;
                    }

                    context.moveTo(p.x, p.y);
                    context.lineTo(p.x + p.vx, p.y + p.vy);
                }

                if (p.dead) {
                    // Remove dead particles, and update the length manually
                    particles.splice(k, 1);
                    particlesLength = particles.length;
                }

            } // End for all particles

            context.stroke();

        } // End if drawParticles

        // lastTime is now
        lastTime = Date.now();

        //var sum = fs.d.reduce(function(a, b) { return a + b; }, 0);
        updateDrawing();
    }

})(window, document);
