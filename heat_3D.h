#ifndef HEAT_3D
#define HEAT_3D

#include <stdbool.h>
#include <stdio.h>
#include "utilities.h"

enum nummethod { FTCS, BE, CN, BEj, BEgs, BEsor };

// 3D discretization
typedef struct {
  double *xx, *yy, *zz;
  int X, Y, Z; // xx, yy, zz are addressed from 1 to X, Y, Z, respectively
} d3D;

// 3D solver prefernces
typedef struct {
  int nx, ny, nz;
  int nsteps, sample, pause;
  double boundary, noise;
  bool periodic, quiet, multigrid;
  enum nummethod method;
  FILE *os, *op;
  double (*source)(long i, long j, long k, const d3D *disc);
  double w; // omega for SOR
  double LX, LY, LZ;
  double alpha, dt;
  double (*init)(double x, double y, double z);
} prefs3D;

void free_d3D(d3D *disc);
d3D *create_d3D(long X, long Y, long Z);

void mgsolve(const prefs3D *p);
void solve(const prefs3D *p);
void ftcs(const prefs3D *p, t3D *d, t3D *s, double Cx, double Cy, double Cz, const d3D *disc);
void cn(const prefs3D *p, t3D *d, t3D *s, double **A, double Cx, double Cy, double Cz, const d3D *disc);
void becs(const prefs3D *p, t3D *d, t3D *s, double **A, const d3D *disc);
void bej(const prefs3D *p, t3D *d, t3D *s, double Cx, double Cy, double Cz, const d3D *disc);
void begs(const prefs3D *p, t3D *d, t3D *s, double Cx, double Cy, double Cz, const d3D *disc);
void besor(const prefs3D *p, t3D *d, t3D *s, double Cx, double Cy, double Cz, const d3D *disc);
void populate_becs_matrix(const prefs3D *p, double **A, long X, long Y, long Z, double Cx, double Cy, double Cz);
void set_constant_boundary(const prefs3D *p, t3D *t);
void set_initial_with_noise(const prefs3D *p, t3D *t, const const d3D *disc, double(*init)(double,double,double));
double gauss3(double x, double y, double z);
double plane_source(long i, long j, long k, const d3D *disc);
void rstrct(double ***d, double ***src, long nh);
void slvsml(const prefs3D *p, double ***solution, double ***rhs);
void fill0(double ***m, long n);
void interp(double ***dst, double ***src, long nh);
void copy(double ***dst, double ***src, long nh);
void relax(const prefs3D *p, double ***sol, double ***rhs, int nf);
void resid(const prefs3D *p, double ***res, double ***sol, double ***rhs, int nf);
void addint(double ***solcur, double ***solprev, double ***res, int nf);

#endif
