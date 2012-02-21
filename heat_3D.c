#include "heat_3D.h"
#include "gauss_elim.h" // gaussian_elimination
#include "utilities.h" // screen, etc.
#include "nrutil.h" // dmatrix, dvector, d3tensor, etc.
#include <signal.h> // signal
#include <stdio.h>  // printf
#include <math.h>   // exp, pow
#define __USE_BSD
#include <unistd.h> // usleep
#include <stdlib.h> // exit

const char *methodnames[] = {
 [FTCS] = "Forward-Time Center-Space",
   [BE] = "Backward Euler",
   [CN] = "Crank-Nicolson",
};

void reset(int s)
{
  printf("\033[m");
  printf("\033[?25h");
  if (s != SIGSTOP)
    exit(EXIT_FAILURE);
}

void solve(const prefs3D *p,
           double Cx, double Cy, double Cz,
           double(*init)(double,double,double))
{
  signal(SIGFPE, reset);  // works
  signal(SIGINT, reset);  // works
  signal(SIGPIPE, reset); // doesn't reset cursor
  screen("\033[2J");        // clear screen
  screen("\033[?25l");      // hide cursor

  double *x, *y, *z, ***T, ***Tnew;
  double **A, **constA;   // in case we use BE or CN
  int Xdim = p->nx+2;
  int Ydim = p->ny+2;
  int Zdim = p->nz+2;
  x = dvector(1, Xdim);
  y = dvector(1, Ydim);
  z = dvector(1, Zdim);
  for (int i = 1; i <= Xdim; i++) x[i] = (i-1) / (double) (p->nx+1);
  for (int j = 1; j <= Ydim; j++) y[j] = (j-1) / (double) (p->ny+1);
  for (int k = 1; k <= Zdim; k++) z[k] = (k-1) / (double) (p->nz+1);

  T      = d3tensor(1, Xdim, 1, Ydim, 1, Zdim);
  Tnew   = d3tensor(1, Xdim, 1, Ydim, 1, Zdim);

  if (!p->periodic) set_constant_boundary(p, T, 1, Xdim, 1, Ydim, 1, Zdim);
  set_initial_with_noise(p, T, 1, Xdim, 1, Ydim, 1, Zdim, x, y, z, init);

  show_d3tensor("T", T, 1, Xdim, 1, Ydim, 1, Zdim);
  usleep(p->pause*2);

  if (p->method == BE || p->method == CN)
  {
    A = dmatrix(1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
    constA = dmatrix(1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
    if (p->method == CN) { Cx /= 2; Cy /= 2; Cz /= 2; }
    populate_becs_matrix(p, constA, Xdim, Ydim, Zdim, Cx, Cy, Cz); // build it once
  }

  for (int n=1; n <= p->nsteps; n++)        // repeat the loop nsteps times
  {
    if (p->method == FTCS)
    {
      ftcs(p, Tnew, T, 1, Xdim, 1, Ydim, 1, Zdim, Cx, Cy, Cz);
    } else if (p->method == BE)
    {
      copy_dmatrix(A, constA, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim); // copy it every time :(
      becs(Tnew, T, 1, Xdim, 1, Ydim, 1, Zdim, A);
    }
    else if (p->method == CN)
    {
      copy_dmatrix(A, constA, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim); // copy it every time :(
      cn(p, Tnew, T, 1, Xdim, 1, Ydim, 1, Zdim, A, Cx, Cy, Cz);
    }

    copy_d3tensor(T, Tnew, 1, Xdim, 1, Ydim, 1, Zdim); // update T to the new values

    if (n%p->sample == 0)
    {
      printf("%s %d\n", methodnames[p->method], n);
      show_d3tensor("T", T, 1, Xdim, 1, Ydim, 1, Zdim);
      usleep(p->pause);
    }
  }

  free_dvector(x, 1, Xdim);
  free_dvector(y, 1, Ydim);
  free_dvector(z, 1, Zdim);
  if (p->method == BE || p->method == CN)
  {
    free_dmatrix(A, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
    free_dmatrix(constA, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
  }
  free_d3tensor(T, 1, Xdim, 1, Ydim, 1, Zdim);
  free_d3tensor(Tnew, 1, Xdim, 1, Ydim, 1, Zdim);
  screen("\033[?25h"); // show cursor$
}

void populate_becs_matrix(const prefs3D *p, double **A, long Xdim, long Ydim, long Zdim, 
                          double Cx, double Cy, double Cz)
{
  // Prepare matrix A from the the Backward Euler discretization using the Cx, Cy, Cz values

  // malloc doesn't guarentee zeroed memory so we zero it out first
  for (int i = 1; i <= Xdim*Ydim*Zdim; i++)
    for (int j = 1; j <= Xdim*Ydim*Zdim; j++)
      A[i][j] = 0;
  if (p->periodic)
  {
    for (int m = 1; m <= Xdim*Ydim*Zdim; m++)
    {
      // Fill main diagonal
      A[m][m] = 2*(Cx+Cy+Cz)+1;
      int left, right;
  
      // Fill off diagonals
      left  = ((m-1)%Zdim == 0)      ? m+Zdim-1 : m-1;
      right = ((m-1)%Zdim == Zdim-1) ? m-Zdim+1 : m+1;
      A[m][left] = A[m][right] = -Cz;
  
      // Fill near bands
      left  = (((m-1)/Zdim)%Ydim == 0)      ? m+Zdim*(Ydim-1) : m-Zdim;
      right = (((m-1)/Zdim)%Ydim == Ydim-1) ? m-Zdim*(Ydim-1) : m+Zdim;
      A[m][left] = A[m][right] = -Cy;
  
      // Fill far bands
      left  = (((m-1)/(Zdim*Ydim)%Xdim) == 0)      ? m+Zdim*Ydim*(Xdim-1) : m-Zdim*Ydim;
      right = (((m-1)/(Zdim*Ydim)%Xdim) == Xdim-1) ? m-Zdim*Ydim*(Xdim-1) : m+Zdim*Ydim;
      A[m][left] = A[m][right] = -Cx;
    }
  } else { // constant boundary
    for (int m = 1; m <= Xdim*Ydim*Zdim; m++)
    {
      // if m is on any boundary, set diagonal to 1 and continue
      if ((m-1)%Zdim == 0 || (m-1)%Zdim == Zdim-1 ||
          ((m-1)/Zdim)%Ydim == 0 || ((m-1)/Zdim)%Ydim == Ydim-1 ||
          ((m-1)/(Zdim*Ydim)%Xdim) == 0 || ((m-1)/(Zdim*Ydim)%Xdim) == Xdim-1)
      {
        A[m][m] = 1;
        continue;
      }
      // Fill main diagonal
      A[m][m] = 2*(Cx+Cy+Cz)+1;
      // Fill off diagonals
      A[m][m-1] = A[m][m+1] = -Cz;
      // Fill near bands
      A[m][m-Zdim] = A[m][m+Zdim] = -Cy;
      // Fill far bands
      A[m][m-Zdim*Ydim] = A[m][m+Zdim*Ydim] = -Cx;
    }
  }
  // show_dmatrix("A", A, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
}

void cn(const prefs3D *p, double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A, double Cx, double Cy, double Cz)
{
  long X = nrh-nrl+1; // include boundaries
  long Y = nch-ncl+1;
  long Z = ndh-ndl+1;
  double *x = dvector(1, X*Y*Z);
  double *b = dvector(1, X*Y*Z);
  // flatten src into b
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
      {
        double *B = &b[1+((i-nrl)*Y+(j-ncl))*Z+(k-ndl)];
        *B = src[i][j][k];
        if (p->periodic) {
          long left, right;
          left = (i == nrl) ? nrh : i-1; right = (i == nrh) ? nrl : i+1;
          *B += Cx*(src[left][j][k] + src[right][j][k] - 2*src[i][j][k]);
          left = (j == ncl) ? nch : j-1; right = (j == nch) ? ncl : j+1;
          *B += Cy*(src[i][left][k] + src[i][right][k] - 2*src[i][j][k]);
          left = (k == ndl) ? ndh : k-1; right = (k == ndh) ? ndl : k+1;
          *B += Cz*(src[i][j][left] + src[i][j][right] - 2*src[i][j][k]);
        } else if (!(i==nrl || i==nrh || j==ncl ||j==nch || k==ndl || k==ndh)) {
          *B += - 2*Cx*src[i][j][k] - 2*Cy*src[i][j][k] - 2*Cz*src[i][j][k];
          *B += (i==nrl) ? 0 : Cx*src[i-1][j][k];
          *B += (i==nrh) ? 0 : Cx*src[i+1][j][k];
          *B += (j==ncl) ? 0 : Cy*src[i][j-1][k];
          *B += (j==nch) ? 0 : Cy*src[i][j+1][k];
          *B += (k==ndl) ? 0 : Cz*src[i][j][k-1];
          *B += (k==ndh) ? 0 : Cz*src[i][j][k+1];
        }
      }
  // solve Ax=b for x, using elimination
  gaussian_elimination(A, x, b, X*Y*Z);
  
  // unflatten x into dst 
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
        dst[i][j][k] = x[1+((i-nrl)*Y+(j-ncl))*Z+(k-ndl)];
  free_dvector(x, 1, X*Y*Z);
  free_dvector(b, 1, X*Y*Z);
}

void becs(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A)
{
  long X = nrh-nrl+1; // include boundaries
  long Y = nch-ncl+1;
  long Z = ndh-ndl+1;
  double *x = dvector(1, X*Y*Z);
  double *b = dvector(1, X*Y*Z);
  // flatten src into b
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
        b[1+((i-nrl)*Y+(j-ncl))*Z+(k-ndl)] = src[i][j][k];

  // solve Ax=b for x, using elimination
  gaussian_elimination(A, x, b, X*Y*Z);

  // unflatten x into dst 
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
        dst[i][j][k] = x[1+((i-nrl)*Y+(j-ncl))*Z+(k-ndl)];

  free_dvector(x, 1, X*Y*Z);
  free_dvector(b, 1, X*Y*Z);
}

void ftcs(const prefs3D *p, double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double Cx, double Cy, double Cz)
{
  for (int i = nrl; i <= nrh; i++)            // calculate new values for all cells
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
      {
        dst[i][j][k] = src[i][j][k];
        if (p->periodic) {
          long left, right;
          left = (i == nrl) ? nrh : i-1; right = (i == nrh) ? nrl : i+1;
          dst[i][j][k] += Cx*(src[left][j][k] + src[right][j][k] - 2*src[i][j][k]);
          left = (j == ncl) ? nch : j-1; right = (j == nch) ? ncl : j+1;
          dst[i][j][k] += Cy*(src[i][left][k] + src[i][right][k] - 2*src[i][j][k]);
          left = (k == ndl) ? ndh : k-1; right = (k == ndh) ? ndl : k+1;
          dst[i][j][k] += Cz*(src[i][j][left] + src[i][j][right] - 2*src[i][j][k]);
        } else if (!(i==nrl || i==nrh || j==ncl ||j==nch || k==ndl || k==ndh)) {
            dst[i][j][k] += Cx*(src[i-1][j][k] + src[i+1][j][k] - 2*src[i][j][k]) +
                            Cy*(src[i][j-1][k] + src[i][j+1][k] - 2*src[i][j][k]) +
                            Cz*(src[i][j][k-1] + src[i][j][k+1] - 2*src[i][j][k]);
        }
      }
}

void set_constant_boundary(const prefs3D *p, double ***T,
                           long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  // first row
  for (int j = ncl; j <= nch; j++)
    for (int k = ndl; k <= ndh; k++)
      T[nrl][j][k] = p->boundary;
  // middle rows
  for (int i = nrl+1; i < nrh; i++)
  {
    // first column, then middle columns, than last column
    for (int k = ndl; k <= ndh; k++) T[i][ncl][k] = p->boundary;
    for (int j = ncl+1; j < nch; j++) T[i][j][ndl] = T[i][j][ndh] = p->boundary;
    for (int k = ndl; k <= ndh; k++) T[i][nch][k] = p->boundary;
  }
  // last row
  for (int j = ncl; j <= nch; j++)
    for (int k = ndl; k <= ndh; k++)
      T[nrh][j][k] = p->boundary;
}

void set_initial_with_noise(const prefs3D *p, double ***T,
                            long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
                            double *x, double *y, double *z,
                            double(*init)(double,double,double))
{
  // if the initial number is I, the result will be between (1-noise)*I and (1+noise)*I
  double A = 1 - p->noise;
  double B = 2*p->noise/RAND_MAX;
  if (!p->periodic) { nrl++; nrh--; ncl++; nch--; ndl++; ndh--; }
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
        T[i][j][k] = init(x[i], y[j], z[k])*(A + B*rand());
}

double gauss3(double x, double y, double z)
{
  return exp(-pow((5*x-2.5), 2))*exp(-pow((5*y-2.5), 2))*exp(-pow((5*z-2.5), 2));
}
