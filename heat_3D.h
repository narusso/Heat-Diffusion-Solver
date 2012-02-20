#ifndef HEAT_3D
#define HEAT_3D

#include <stdbool.h>

enum method { FTCS, BE, CN };

void discretize(int nx, int ny, int nz, int nsteps, int sample, int pause,
                double Cx, double Cy, double Cz,
                double(*f)(double,double,double), double max, double boundary, bool periodic, enum method);
void ftcs(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double Cx, double Cy, double Cz, bool periodic);
void cn(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A, double Cx, double Cy, double Cz, bool periodic);
void becs(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A, bool periodic);
void populate_becs_matrix(double **A, long Xdim, long Ydim, long Zdim,
                          double Cx, double Cy, double Cz, bool periodic);
void set_constant_boundary(double ***T,
                           long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
                           double value);
void set_initial_with_noise(double ***T,
                            long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
                            double *x, double *y, double *z,
                            double(*f)(double,double,double), double max);
double gauss3(double x, double y, double z);

#endif
