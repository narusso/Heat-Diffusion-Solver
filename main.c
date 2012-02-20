// Forward-Time, Centered Space, and Backward Euler, 3D Heat Diffusion in 1m^3 of copper
#include "heat_3D.h"
#include <stdio.h>  // fprintf
#include <unistd.h> // getpid
#include <time.h>   // time
#include <stdlib.h> // exit
#include <getopt.h> // getopt
#include <string.h> // strlen
#include <stdbool.h> // bool
#include <assert.h> // assert
void usage(char *name);

int main(int argc, char *argv[])
{
  double LX, LY, LZ, alpha, dx, dy, dz, dt, Cx, Cy, Cz, noise, boundary;
  int nx, ny, nz, nsteps, sample, pause;
  enum method meth = BE;
  bool periodic;

  // Default values
  LX = LY = LZ = 1;                      // length of 1D object in m along each axis
  nx = ny = nz = 3;                      // number of divisions along each axis
  nsteps = 6;                            // number of time steps
  sample = 2;                            // how often to show output
  pause  = 10000;                        // how long to pause after showing output
  alpha  = 1.1234e-4;                    // diffusivity of copper in m^2/s
  dt     = .003;                         // length of one time step in seconds
  noise  = 0.2;                          // max noise to add/subtract from initial values
  boundary = 0;                          // constant boundary condition
  periodic = false;                      // periodic boundary condition

  int opt;
  while ((opt = getopt(argc, argv, "X:Y:Z:x:y:z:n:s:p:a:t:r:b:m:")) != -1)
  {
    switch (opt) {
      case 'X': LX = atof(optarg); break;
      case 'Y': LY = atof(optarg); break;
      case 'Z': LZ = atof(optarg); break;
      case 'x': nx = atoi(optarg); break;
      case 'y': ny = atoi(optarg); break;
      case 'z': nz = atoi(optarg); break;
      case 'n': nsteps = atoi(optarg); break;
      case 's': sample = atoi(optarg); break;
      case 'p': pause = atof(optarg)*1000*1000; break;
      case 'a': alpha = atof(optarg); break;
      case 't': dt = atof(optarg); break;
      case 'r': noise = atof(optarg); break;
      case 'b':
        if (strcmp(optarg, "p") == 0) periodic = true;
        else boundary = atof(optarg);
        break;
      case 'm':
        if (strcmp(optarg, "FTCS") == 0) { meth = FTCS; break; }
        else if (strcmp(optarg, "BE") == 0) { meth = BE; break; }
        else if (strcmp(optarg, "CN") == 0) { meth = CN; break; }
      default:
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  assert(LX > 0); assert(LY > 0); assert(LZ > 0);
  assert(nx > 0); assert(ny > 0); assert(nz > 0);
  assert(sample > 0); assert(alpha > 0);
  assert(dt > 0); 

  dx = LX/nx; dy = LY/ny; dz = LZ/nz;    // length of one segment in m along z axis
  Cx = alpha*dt/(dx*dx);
  Cy = alpha*dt/(dy*dy);
  Cz = alpha*dt/(dz*dz);

  srand(getpid()*time(NULL));
  solve(nx, ny, nz, nsteps, sample, pause, Cx, Cy, Cz, gauss3, noise, boundary, periodic, meth);
  exit(EXIT_SUCCESS);
}

void usage(char *name)
{
  const char *option[] = {
    "[-X length in x dimension (in meters) ]", "[-Y length in y dimension (in meters) ]",
    "[-Z length in z dimension (in meters) ]", "[-x number of divisions along x dimension ]",
    "[-y number of divisions along y dimension ]", "[-z number of divisions along z dimension ]",
    "[-n number of time steps to calculate ]", "[-s how many time steps between reports ]",
    "[-p how long to pause reports (in seconds) ]", "[-a diffusivity constant (in m/s^2) ]",
    "[-t length of time step (in seconds) ]", "[-r ratio of noise applied to initial condition (0=none) ]",
    "[-b value of constant boundary condition (or p for periodic) ]", "[-m method to use (FTCS, BE, or CN) ]"
  };
  int indentation = strlen("Usage:  ") + strlen(name);
  fprintf(stderr, "Usage: %s ", name);
  fprintf(stderr, "%s\n", option[0]);
  for (int i = 1; i < (sizeof(option)/sizeof(option[0])); i++)
  {
    for (int j = 0; j < indentation; j++) fprintf(stderr, " ");
    fprintf(stderr, "%s\n", option[i]);
  }
}
