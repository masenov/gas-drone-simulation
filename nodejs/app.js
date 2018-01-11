var express = require('express');
var app = express();
var http = require('http').Server(app);
var io = require('socket.io')(http);
var fsl = require('fs');
eval(fsl.readFileSync('fluidsolver_ss.js')+'');

app.use(express.static(__dirname + "/"));
app.get('/', function(req, res){
    res.sendFile(__dirname + '/index.html');
});

http.listen(3000, function(){
    console.log('listening on *:3000');
});


var NUM_OF_CELLS = 256, // Number of cells (not including the boundary)
    VIEW_SIZE = 640,    // View size (square)
    FPS = 60;           // Frames per second
var CELL_SIZE = VIEW_SIZE / NUM_OF_CELLS,  // Size of each cell in pixels
    CELL_SIZE_CEIL = Math.ceil(CELL_SIZE); // Size of each cell in pixels (ceiling)

/**
 * A simple particle class.
 *
 * @param x {Number}
 * @param y {Number}
 * @constructor
 */
function Particle(x, y) {
    this.x = x;
    this.y = y;
    this.vx = 0;
    this.vy = 0;
    this.age = 0;
    this.dead = false;
}

Particle.TIME_TO_LIVE = 5; // Time to live in seconds


// Create the fluid solver
var fs = new FluidSolver(NUM_OF_CELLS);
fs.resetVelocity();

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


function toRadians (angle) {
    return angle * (Math.PI / 180);
}

/**
 * @param min {Number}
 * @param max {Number}
 * @returns {Number}
 */
function getRandom(min, max) {
    return min + Math.random() * ((max + 1) - min);
}




    /**
     * Update loop
     */
    function update(/*time*/) {
        var i, j, k, l, m, dx, dy, color, cellIndex, deltaTime,
            r, g, b, u, v, density, pxIdx, pxX, pxY,
            invMaxColor = 1.0 / 255;

        fs.dOld[fs.I(fs.gasLocationX, fs.gasLocationY)] = fs.gasRelease;
        console.log(fs.d[fs.I(10,10)]);
        deltaTime = (Date.now() - lastTime) / 1000;

        // Step the fluid simulation
        fs.velocityStep();
        fs.densityStep();

        var sum = fs.d.reduce(function(a, b) { return a + b; }, 0);
        console.log(sum);
        // Artificial wind field
        var du = Math.sin(toRadians(fs.windDirection-180))*fs.windSpeed,
            dv = Math.cos(toRadians(fs.windDirection-180))*fs.windSpeed;
        var ii,J;
        // Add the mouse velocity to cells above, below, to the left, and to the right as well.
        var np = fs.windLocations;
        var acc = NUM_OF_CELLS/np;
        for (i = 0; i<np; i++) {
            for (j = 0; j < np; j++) {
                ii = i*acc;
                jj = j*acc;
                fs.uOld[fs.I(ii, jj)] = du;
                fs.vOld[fs.I(ii, jj)] = dv;

            }
        }
        return fs;

    } // End update()

io.on('connection', function(socket){
    socket.on('update simulation', function(msg){
        console.log("Received request!");
        newfs = update();
        socket.emit('receive simulation data', newfs);
    });
});

io.on('connection', function(socket){
    socket.on('chat message', function(msg){
        console.log('message: ' + msg);
        io.emit('chat message', 'test');
    });
});

