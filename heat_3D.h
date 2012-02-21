#ifndef HEAT_3D
#define HEAT_3D

#include <stdbool.h>

enum nummethod { FTCS, BE, CN };

// 3D Temp array and bounds
typedef struct {
  double ***T;
  long nrl, nrh, ncl, nch, ndl, ndh;
} temp3D;

// 3D solver prefernces
typedef struct {
  int nx, ny, nz;
  int nsteps, sample, pause;
  double boundary, noise;
  bool periodic;
  enum nummethod method;
} prefs3D;


void solve(const prefs3D *p, int nx, int ny, int nz, int nsteps, int sample, int pause,
           double Cx, double Cy, double Cz,
           double(*f)(double,double,double));
void ftcs(const prefs3D *p, double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double Cx, double Cy, double Cz);
void cn(const prefs3D *p, double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A, double Cx, double Cy, double Cz);
void becs(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A);
void populate_becs_matrix(const prefs3D *p, double **A, long Xdim, long Ydim, long Zdim,
                          double Cx, double Cy, double Cz);
void set_constant_boundary(const prefs3D *p, double ***T,
                           long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void set_initial_with_noise(const prefs3D *p, double ***T,
                            long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
                            double *x, double *y, double *z,
                            double(*f)(double,double,double));
double gauss3(double x, double y, double z);

#endif
