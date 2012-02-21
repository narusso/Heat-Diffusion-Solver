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
  // parameters
  prefs3D ps, *p; p = &ps;
  // parameters not directly used by solver
  double LX, LY, LZ, alpha, dt;

  // default parameter values
  p->nx = p->ny = p->nz = 3;             // number of divisions along each axis
  p->nsteps = 6;                         // number of time steps
  p->sample = 2;                         // how often to show output
  p->pause  = 10000;                     // how long to pause after showing output
  p->noise  = 0.2;                       // max noise to add/subtract from initial values
  p->boundary = 0;                       // constant boundary condition
  p->periodic = false;                   // periodic boundary condition
  p->method = BE;
  LX = LY = LZ = 1;                      // length of 1D object in m along each axis
  alpha  = 1.1234e-4;                    // diffusivity of copper in m^2/s
  dt     = .003;                         // length of one time step in seconds

  int opt;
  while ((opt = getopt(argc, argv, "X:Y:Z:x:y:z:n:s:p:a:t:r:b:m:")) != -1)
  {
    switch (opt) {
      case 'X': LX = atof(optarg); break;
      case 'Y': LY = atof(optarg); break;
      case 'Z': LZ = atof(optarg); break;
      case 'x': p->nx = atoi(optarg); break;
      case 'y': p->ny = atoi(optarg); break;
      case 'z': p->nz = atoi(optarg); break;
      case 'n': p->nsteps = atoi(optarg); break;
      case 's': p->sample = atoi(optarg); break;
      case 'p': p->pause = atof(optarg)*1000*1000; break;
      case 'a': alpha = atof(optarg); break;
      case 't': dt = atof(optarg); break;
      case 'r': p->noise = atof(optarg); break;
      case 'b':
        if (strcmp(optarg, "p") == 0) p->periodic = true;
        else p->boundary = atof(optarg);
        break;
      case 'm':
        if (strcmp(optarg, "FTCS") == 0) { p->method = FTCS; break; }
        else if (strcmp(optarg, "BE") == 0) { p->method = BE; break; }
        else if (strcmp(optarg, "CN") == 0) { p->method = CN; break; }
      default:
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  assert(LX > 0); assert(LY > 0); assert(LZ > 0);
  assert(p->nx > 0); assert(p->ny > 0); assert(p->nz > 0);
  assert(p->sample > 0); assert(alpha > 0); assert(p->pause > 0);
  assert(dt > 0); 

  // derived constants
  double dx, dy, dz, Cx, Cy, Cz;
  dx = LX/p->nx; dy = LY/p->ny; dz = LZ/p->nz;
  Cx = alpha*dt/(dx*dx);
  Cy = alpha*dt/(dy*dy);
  Cz = alpha*dt/(dz*dz);

  srand(getpid()*time(NULL));
  solve(p, Cx, Cy, Cz, gauss3);
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
