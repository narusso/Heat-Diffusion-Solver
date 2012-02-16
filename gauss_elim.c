#include "gauss_elim.h"

void gaussian_elimination(double **A, double *x, double *b, int n)
{
  upper_triangulate(A,b,n);
  back_sub(A,x,b,n);
}

void swap_rows(double **A, double *b, int j, int best, int n)
{
  double *r, t;

  r = (double *)malloc((size_t) (n*sizeof(double)));
  if (!r) nrerror("allocation failure in swap_rows()");

  memcpy(r,           &A[j][1],    n*sizeof(double));
  memcpy(&A[j][1],    &A[best][1], n*sizeof(double));
  memcpy(&A[best][1], r,           n*sizeof(double));

  t = b[j];
  b[j] = b[best];
  b[best] = t;

  free(r);
}

void upper_triangulate(double **A, double *b, int n)
{
  double r;                 /* obviously, this is the scale factor */
  int i,j,k;
  for (j=1;j<n;++j){       /* loop over first n-1 columns */
    { // exchange row to maximize pivot
      int h,best;
      double pivot;
      best = j; pivot = A[j][j];
      for(h=j+1;h<=n;++h){
        if (fabs(A[h][j]) > fabs(pivot)){
          best = h; pivot = A[h][j];
        }
      }
      if (best != j) swap_rows(A, b, j, best, n);
    }
    for(i=j+1;i<=n;++i){    /* loop over rows */
      if (A[i][j] != 0.0){  /* if element not already zeroed out */
        r = A[i][j]/A[j][j];/* multiple of pivot row */
        for (k=1;k<=n;++k){ /* for row element carry out subtraction */
          A[i][k] = A[i][k] - r*A[j][k];
        }
        b[i] = b[i] - r*b[j]; /* for b too */
      }
    }
  }
}

void back_sub(double **A, double *x, double *b, int m)
{
  int i,j;
  x[m] = b[m]/A[m][m];

  for (i=m-1;i>=1;--i){
    x[i] = b[i];
    for (j=i+1;j<=m;++j){
      x[i] -= A[i][j]*x[j];
    }
    x[i]/=A[i][i];
  }
}
