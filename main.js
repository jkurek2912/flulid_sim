class Sim {
    constructor(numX, numY, h) {
        // boundary cells
        this.numX = numX + 2;
        this.numY = numY + 2;
        this.numCells = this.numX * this.numY;
         
        // cell size
        this.h = h;
        
        // velocity field
        this.u = new Float32Array(this.numCells);
        this.v = new Float32Array(this.numCells);

        // array for advection
        this.newU = new Float32Array(this.numCells);
        this.newV = new Float32Array(this.numCells);

        // array to keep track of state (solid = 0, fluid = 1)
        this.s = new Float32Array(this.numCells);

        // smoke/dye concentration
        this.m = new Float32Array(this.numCells);
        this.newM = new Float32Array(this.numCells);

        this.s.fill(1);
        this.m.fill(0);
    }

    addForces(dt, gravity) {
        var n = this.numY;

        // loop through non-boundary cells
        for (var i = 1; i < this.numX; i++) {
            for (var j = 1; j < this.numY - 1; j++) {
            
                // apply forces to fluid cells that have fluid below them
                // use i*n+j converts 2d coordinates to 1d state array index (row major ordering)
                if (this.s[i*n + j] != 0 && this.s[i*n + j - 1] != 0) {

                    // add gravity to vertical velocity
                    this.v[i*n + j] += gravity *dt;
                }
            }
        }
    }

    correctDivergence(numIterations, dt, overRelaxation) {
        var n = this.numY;

        // iterate multiple times for better approximation
        for (var iter = 0; iter < numIterations; iter++) {

            // loop through non-boundary cells
            for (var i = 1; i < this.numX - 1; i++) {
                for (var j = 1; j < this.numY - 1; j++) {

                    // skip solid cells
                    if (this.s[i*n + j] == 0) {
                        continue;
                    }

                    // check which neighboring cells are fluid
                    var nx1 = this.s[(i+1)*n + j];  // right
                    var nx2 = this.s[(i-1)* n + j]; // left
                    var ny1 = this.s[i*n + j + 1];    // top
                    var ny2 = this.s[i*n + j - 1];    // bottom

                    var s = nx1 + nx2 + ny1 + ny2;

                    // skip cells with not fluid neighbors
                    if (s == 0) {
                        continue;
                    }

                    // calculate divergence
                    var div = this.u[(i+1)*n + j] - this.u[i*n + j] + 
                         this.v[i*n + j+1] - this.v[i*n + j];

                    // calculate pressure
                    var p = -div / s * overRelaxation;

                    // apply presure adjustment to cells
                    this.u[(i+1)*n + j] += nx1 * p;     // Add to right face
                    this.u[i*n + j] -= nx2 * p;         // Subtract from left face
                    this.v[i*n + j+1] += ny1 * p;       // Add to top face
                    this.v[i*n + j] -= ny2 * p;         // Subtract from bottom face
                    
                }
            }
        }
    }

    fillBoundaryVelocities() {
        var n = this.numY;

        // copy horizontal velocities to boundary cells
        for (var i = 0; i < this.numX; i++) {
            
            // bottom boundary
            this.u[i*n] = this.u[i*n + 1];
            
            // top boundary
            this.u[i*n + this.numY - 1] = this.u[i*n + this.numY - 2];
        }

        // copy vertical velocities to boundary cells
        for (var j = 0; j < this.numY; j++) {
            
            // left boundary
            this.v[j] = this.v[n + j];
            
            // right boundary
            this.v[(this.numX - 1)*n + j] = this.v[(this.numX - 2)*n + j];
        }
    }


}