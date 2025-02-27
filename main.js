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

        this.s.fill(1.0);
        this.m.fill(1.0);
    }

    addForces(dt, gravity) {
        var n = this.numY;

        // loop through all grid cells except for boundary
        for (var i = 1; i < this.numX; i++) {
            for (var j = 1; j < this.numY; j++) {
            
                // apply forces to fluid cells that have fluid below them
                // use i*n+j converts 2d coordinates to 1d state array index (row major ordering)
                if (this.s[i*n+j] != 0.0 && this.s[i*n+j-1] != 0.0) {

                    // add gravity to vertical velocity
                    this.v[i*n+j] += gravity *dt;
                }
            }
        }
    }


}