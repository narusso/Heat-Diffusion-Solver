#ifndef GAUSS_ELIM
#define GAUSS_ELIM
#include "nrutil.h"
#include <stdlib.h> // malloc and free
#include <string.h> // memcpy
#include <math.h>   // fabs
#include <assert.h> // assert

void gaussian_elimination(double **A, double *x, double *b, int n);
void upper_triangulate(double **A, double *b, int n);
void back_sub(double **A, double *x, double *b, int n);

#endif
