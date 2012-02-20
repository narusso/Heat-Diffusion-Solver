#ifndef GAUSS_ELIM
#define GAUSS_ELIM

void gaussian_elimination(double **A, double *x, double *b, int n);
void upper_triangulate(double **A, double *b, int n);
void back_sub(double **A, double *x, double *b, int n);

#endif
