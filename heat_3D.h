#ifndef HEAT_3D
#define HEAT_3D

#include <stdbool.h>
#include "utilities.h"

enum nummethod { FTCS, BE, CN };

// 3D solver prefernces
typedef struct {
  int nx, ny, nz;
  int nsteps, sample, pause;
  double boundary, noise;
  bool periodic;
  enum nummethod method;
} prefs3D;

void solve(const prefs3D *p, double Cx, double Cy, double Cz, double(*init)(double,double,double));
void ftcs(const prefs3D *p, temp3D *d, temp3D *s, double Cx, double Cy, double Cz);
void cn(const prefs3D *p, temp3D *d, temp3D *s, double **A, double Cx, double Cy, double Cz);
void becs(temp3D *d, temp3D *s, double **A);
void populate_becs_matrix(const prefs3D *p, double **A, long X, long Y, long Z, double Cx, double Cy, double Cz);
void set_constant_boundary(const prefs3D *p, temp3D *t);
void set_initial_with_noise(const prefs3D *p, temp3D *t, double *x, double *y, double *z, double(*init)(double,double,double));
double gauss3(double x, double y, double z);

#endif
