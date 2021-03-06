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
void show_params(FILE *f, int argc, char **argv);
void usage(char *name);

int main(int argc, char *argv[])
{
  // parameters
  prefs3D ps, *p; p = &ps;
  char *out_soln = NULL;                  // filename for solution data
  char *out_perf = NULL;                  // filename for performance data

  // default parameter values
  p->nx = p->ny = p->nz = 3;              // number of internal points along each axis
                                          // nx+2 is total number of points
                                          // nx+1 is number of intervals
                                          // LX/(nx+1) is size of interval
  p->nsteps    = 6;                       // number of time steps
  p->sample    = 2;                       // how often to show output
  p->pause     = 10000;                   // how long to pause after showing output
  p->noise     = 0.2;                     // max noise to add/subtract from initial values
  p->boundary  = 0;                       // constant boundary condition
  p->periodic  = false;                   // periodic boundary condition
  p->method    = BE;                      // method for solving Ax=b
  p->quiet     = false;                   // whether to print data to stdout
  p->multigrid = false;                   // whether to use multigrid w/BE
  p->os        = NULL;                    // FILE* for solution data
  p->op        = NULL;                    // FILE* for performance data
  p->source    = plane_source;            // only one source function available
  p->LX = p->LY = p->LZ = 1;              // length of 1D object in m along each axis
  p->alpha     = 1.1234e-4;               // diffusivity of copper in m^2/s
  p->dt        = .003;                    // length of one time step in seconds
  p->w         = 1.65;                    // omega for SOR
  p->init      = gauss3;                  // initialize domain with 3D Gaussian

  int opt;
  while ((opt = getopt(argc, argv, "X:Y:Z:x:y:z:n:s:p:a:t:r:b:m:o:O:w:qg")) != -1)
  {
    switch (opt) {
      case 'X': p->LX = atof(optarg); break;
      case 'Y': p->LY = atof(optarg); break;
      case 'Z': p->LZ = atof(optarg); break;
      case 'x': p->nx = atoi(optarg); break;
      case 'y': p->ny = atoi(optarg); break;
      case 'z': p->nz = atoi(optarg); break;
      case 'n': p->nsteps = atoi(optarg); break;
      case 's': p->sample = atoi(optarg); break;
      case 'p': p->pause = atof(optarg)*1000*1000; break;
      case 'a': p->alpha = atof(optarg); break;
      case 't': p->dt = atof(optarg); break;
      case 'r': p->noise = atof(optarg); break;
      case 'w': p->w = atof(optarg); break;
      case 'q': p->quiet = true; break;
      case 'g': p->multigrid = true; break;
      case 'b':
        if (strcmp(optarg, "p") == 0) p->periodic = true;
        else p->boundary = atof(optarg);
        break;
      case 'm':
        if (strcmp(optarg, "FTCS") == 0) { p->method = FTCS; break; }
        else if (strcmp(optarg, "BE") == 0) { p->method = BE; break; }
        else if (strcmp(optarg, "CN") == 0) { p->method = CN; break; }
        else if (strcmp(optarg, "BEj") == 0) { p->method = BEj; break; }
        else if (strcmp(optarg, "BEgs") == 0) { p->method = BEgs; break; }
        else if (strcmp(optarg, "BEsor") == 0) { p->method = BEsor; break; }
        else { usage(argv[0]); exit(EXIT_FAILURE); }
      case 'o':
        if (strcmp(optarg, ".") == 0) out_soln = "/tmp/heat_3D_soln.dat";
        else out_soln = optarg;
        break;
      case 'O':
        if (strcmp(optarg, ".") == 0) out_perf = "/tmp/heat_3D_perf.dat";
        else out_perf = optarg;
        break;
      default:
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  assert(p->LX > 0); assert(p->LY > 0); assert(p->LZ > 0);
  assert(p->nx > 0); assert(p->ny > 0); assert(p->nz > 0);
  assert(p->sample > 0); assert(p->alpha > 0); assert(p->pause >= 0);
  assert(p->dt > 0); 

  // To use the multigrid method, number of intervals must be a power of 2
  if (p->multigrid)
  {
    int n; // number of intervals
    int e; // exponent
    for (n=p->nx+1, e=0; n >>= 1; e++) ;
    if (p->nx+1 != (1L << e)) { fprintf(stderr,"nx+1 must be a power of 2 for multigrid method\n"); exit(1); }
    for (n=p->ny+1, e=0; n >>= 1; e++) ;
    if (p->ny+1 != (1L << e)) { fprintf(stderr,"ny+1 must be a power of 2 for multigrid method\n"); exit(1); }
    for (n=p->nz+1, e=0; n >>= 1; e++) ;
    if (p->nz+1 != (1L << e)) { fprintf(stderr,"nz+1 must be a power of 2 for multigrid method\n"); exit(1); }
    if (p->nx != p->ny || p->ny != p->nz) { fprintf(stderr,"For multigrid method, nx==ny==nz\n"); exit(1); }
  }

  srand(getpid()*time(NULL));
  if (out_soln)
  {
    p->os = fopen(out_soln, "w");
    if (!p->os)
    {
      fprintf(stderr, "Couldn't open %s for writing\n", out_soln);
      exit(EXIT_FAILURE);
    }
    show_params(p->os, argc, argv);
  }
  if (out_perf)
  {
    p->op = fopen(out_perf, "w");
    if (!p->op)
    {
      fprintf(stderr, "Couldn't open %s for writing\n", out_soln);
      exit(EXIT_FAILURE);
    }
    show_params(p->op, argc, argv);
  }
  if (p->multigrid) {
    mgsolve(p);
  } else {
    solve(p);
  }
  if (p->os) fclose(p->os);
  if (p->op) fclose(p->op);
  if (!p->quiet) show_params(stdout, argc, argv);
  exit(EXIT_SUCCESS);
}

void show_params(FILE *f, int argc, char **argv)
{
    // Output parameters into data file for reproducibility
    fprintf(f, "# ");
    for (int i = 0; i < argc; i++) fprintf(f, "%s ", argv[i]);
    fprintf(f, "\n");
}

void usage(char *name)
{
  const char *option[] = {
    "[-X length in x dimension (in meters) ]", "[-Y length in y dimension (in meters) ]",
    "[-Z length in z dimension (in meters) ]", "[-x number of internal points along x dimension ]",
    "[-y number of internal points along y dimension ]", "[-z number of internal points along z dimension ]",
    "[-n number of time steps to calculate ]", "[-s how many time steps between reports ]",
    "[-p how long to pause reports (in seconds) ]", "[-a diffusivity constant (in m/s^2) ]",
    "[-t length of time step (in seconds) ]", "[-r ratio of noise applied to initial condition (0=none) ]",
    "[-b value of constant boundary condition (or p for periodic) ]",
    "[-m method (FTCS BE CN BEj BEgs BEsor) ]",
    "[-o filename for plottable solution data ]", "[-O filename for plottable performance data ]",
    "[-q suppress normal output ]",
    "[-g use multigrid/BE technique ]",
  };
  int indentation = strlen("Usage:  ") + strlen(name);
  fprintf(stderr, "Usage: %s ", name);
  fprintf(stderr, "%s\n", option[0]);
  for (int i = 1; i < (sizeof(option)/sizeof(option[0])); i++)
  {
    for (int j = 0; j < indentation; j++) fprintf(stderr, " ");
    fprintf(stderr, "%s\n", option[i]);
  }
  fprintf(stderr, "Specify . as the filename to use a default value\n");
}
