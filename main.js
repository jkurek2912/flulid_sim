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
    this.m.fill(0.0);
  }

  addForces(dt, gravity) {
    var n = this.numY;

    // loop through non-boundary cells
    for (var i = 1; i < this.numX; i++) {
      for (var j = 1; j < this.numY - 1; j++) {
        // apply forces to fluid cells that have fluid below them
        // use i*n+j converts 2d coordinates to 1d state array index (row major ordering)
        if (this.s[i * n + j] != 0.0 && this.s[i * n + j - 1] != 0.0) {
          // add gravity to vertical velocity
          this.v[i * n + j] += gravity * dt;
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
          if (this.s[i * n + j] == 0) {
            continue;
          }

          // check which neighboring cells are fluid
          var nx1 = this.s[(i + 1) * n + j]; // right
          var nx2 = this.s[(i - 1) * n + j]; // left
          var ny1 = this.s[i * n + j + 1]; // top
          var ny2 = this.s[i * n + j - 1]; // bottom

          var s = nx1 + nx2 + ny1 + ny2;

          // skip cells with not fluid neighbors
          if (s == 0) {
            continue;
          }

          // calculate divergence
          var div =
            this.u[(i + 1) * n + j] -
            this.u[i * n + j] +
            this.v[i * n + j + 1] -
            this.v[i * n + j];

          // calculate pressure
          var p = (-div / s) * overRelaxation;

          // apply presure adjustment to cells
          this.u[(i + 1) * n + j] += nx1 * p; // Add to right face
          this.u[i * n + j] -= nx2 * p; // Subtract from left face
          this.v[i * n + j + 1] += ny1 * p; // Add to top face
          this.v[i * n + j] -= ny2 * p; // Subtract from bottom face
        }
      }
    }
  }

  fillBoundaryVelocities() {
    var n = this.numY;

    // copy horizontal velocities to boundary cells
    for (var i = 0; i < this.numX; i++) {
      // bottom boundary
      this.u[i * n] = this.u[i * n + 1];

      // top boundary
      this.u[i * n + this.numY - 1] = this.u[i * n + this.numY - 2];
    }

    // copy vertical velocities to boundary cells
    for (var j = 0; j < this.numY; j++) {
      // left boundary
      this.v[0 * n + j] = this.v[1 * n + j];

      // right boundary
      this.v[(this.numX - 1) * n + j] = this.v[(this.numX - 2) * n + j];
    }
  }

  sampleGridField(x, y, field) {
    var n = this.numY;
    var h = this.h;
    var h1 = 1.0 / h;
    var h2 = 0.5 * h;

    // set x and y into valid range
    x = Math.max(Math.min(x, this.numX * h), h);
    y = Math.max(Math.min(y, this.numY * h), h);

    var dx = 0.0;
    var dy = 0.0;

    // adjust sampling based on field type
    var f;
    if (field == 0) {
      // horizontal velocity
      f = this.u;
      dy = h2;
    } else if (field == 1) {
      // vertical velocity
      f = this.v;
      dx = h2;
    } else if (field == 2) {
      // smoke/density
      f = this.m;
      dx = h2;
      dy = h2;
    }

    // calculate grid indices and interpolate weights
    var x0 = Math.min(Math.floor((x - dx) * h1), this.numX - 1);
    var tx = (x - dx - x0 * h) * h1;
    var x1 = Math.min(x0 + 1, this.numX - 1);

    var y0 = Math.min(Math.floor((y - dy) * h1), this.numY - 1);
    var ty = (y - dy - y0 * h) * h1;
    var y1 = Math.min(y0 + 1, this.numY - 1);

    // bilinear interpolation weights
    var sx = 1.0 - tx;
    var sy = 1.0 - ty;

    var val =
      sx * sy * f[x0 * n + y0] +
      tx * sy * f[x1 * n + y0] +
      tx * ty * f[x1 * n + y1] +
      sx * ty * f[x0 * n + y1];

    return val;
  }

  advectVelocity(dt) {
    // save current velocity field
    this.newU.set(this.u);
    this.newV.set(this.v);

    var n = this.numY;
    var h = this.h;
    var h2 = 0.5 * h;

    // field constants for sampleGridField
    var U_FIELD = 0;
    var V_FIELD = 1;

    // adjust velocity fields
    for (var i = 1; i < this.numX; i++) {
      for (var j = 1; j < this.numY; j++) {
        // advect horizontal velocity
        if (
          this.s[i * n + j] != 0.0 &&
          this.s[(i - 1) * n + j] != 0.0 &&
          j < this.numY - 1
        ) {
          var x = i * h;
          var y = j * h + h2;

          // get velocity field at this position
          var u = this.u[i * n + j];
          var v = this.avgV(i, j);

          // backtrace and clamp
          x = x - dt * u;
          y = y - dt * v;
          x = Math.max(Math.min(x, this.numX * h), h);
          y = Math.max(Math.min(y, this.numY * h), h);

          // sample velocity field from previous position
          u = this.sampleGridField(x, y, U_FIELD);
          this.newU[i * n + j] = u;
        }

        // advect vertical velocity
        if (
          this.s[i * n + j] != 0.0 &&
          this.s[i * n + j - 1] != 0.0 &&
          i < this.numX - 1
        ) {
          var x = i * h + h2;
          var y = j * h;

          // get velocity field at this position
          var u = this.avgU(i, j);
          var v = this.v[i * n + j];

          // backtrace and clamp
          x = x - dt * u;
          y = y - dt * v;
          x = Math.max(Math.min(x, this.numX * h), h);
          y = Math.max(Math.min(y, this.numY * h), h);

          v = this.sampleGridField(x, y, V_FIELD);
          this.newV[i * n + j] = v;
        }
      }
    }
    // update velocity fields
    this.u.set(this.newU);
    this.v.set(this.newV);
  }

  avgU(i, j) {
    var n = this.numY;
    var u =
      (this.u[i * n + j - 1] +
        this.u[i * n + j] +
        this.u[(i + 1) * n + j - 1] +
        this.u[(i + 1) * n + j]) *
      0.25;
    return u;
  }

  avgV(i, j) {
    var n = this.numY;
    var v =
      (this.v[(i - 1) * n + j] +
        this.v[i * n + j] +
        this.v[(i - 1) * n + j + 1] +
        this.v[i * n + j + 1]) *
      0.25;
    return v;
  }

  advectSmoke(dt) {
    // save current smoke field
    this.newM.set(this.m);

    var n = this.numY;
    var h = this.h;
    var h2 = 0.5 * h;

    // define field constant for sampleGridField
    var S_FIELD = 2;

    // advect smoke field
    for (var i = 1; i < this.numX; i++) {
      for (var j = 1; j < this.numY; j++) {
        if (this.s[i * n + j] != 0.0) {
          // get velocity at cell center
          var u = (this.u[i * n + j] + this.u[(i + 1) * n + j]) * 0.5;
          var v = (this.v[i * n + j] + this.v[i * n + j + 1]) * 0.5;

          // trace back
          var x = i*h + h2 - dt*u;
          var y = j*h + h2 - dt*v;

          // sample smoke from previous position
          this.newM[i*n + j] = this.sampleGridField(x, y, S_FIELD);
        }
      }
    }

    // update smoke field
    this.m.set(this.newM);
  }
}
