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

const char *methods[] = {
 [FTCS] = "Forward-Time Center-Space",
   [BE] = "Backward Euler",
   [CN] = "Crank-Nicolson",
};

void reset(int s)
{
  printf("[m");
  printf("[?25h");
  if (s != SIGSTOP)
    exit(EXIT_FAILURE);
}

void discretize(int nx, int ny, int nz, int nsteps, int sample, int pause,
                double Cx, double Cy, double Cz,
                double(*init)(double,double,double), double max, double boundary, bool periodic, enum method m)
{
  signal(SIGFPE, reset);  // works
  signal(SIGINT, reset);  // works
  signal(SIGPIPE, reset); // doesn't reset cursor
  screen("[2J");        // clear screen
  screen("[?25l");      // hide cursor

  double *x, *y, *z, ***T, ***Tnew;
  double **A, **constA;   // in case we use BE or CN
  int Xdim = nx+2;
  int Ydim = ny+2;
  int Zdim = nz+2;
  x = dvector(1, Xdim);
  y = dvector(1, Ydim);
  z = dvector(1, Zdim);
  for (int i = 1; i <= Xdim; i++) x[i] = (i-1) / (double) (nx+1);
  for (int j = 1; j <= Ydim; j++) y[j] = (j-1) / (double) (ny+1);
  for (int k = 1; k <= Zdim; k++) z[k] = (k-1) / (double) (nz+1);

  T      = d3tensor(1, Xdim, 1, Ydim, 1, Zdim);
  Tnew   = d3tensor(1, Xdim, 1, Ydim, 1, Zdim);

  set_initial_with_noise(T, 1, Xdim, 1, Ydim, 1, Zdim, x, y, z, init, max);
  if (!periodic) set_constant_boundary(T, 1, Xdim, 1, Ydim, 1, Zdim, boundary);

  show_d3tensor("T", T, 1, Xdim, 1, Ydim, 1, Zdim);
  usleep(pause*2);

  if (m == BE || m == CN)
  {
    A = dmatrix(1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
    constA = dmatrix(1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
    if (m == CN) { Cx /= 2; Cy /= 2; Cz /= 2; }
    populate_becs_matrix(constA, Xdim, Ydim, Zdim, Cx, Cy, Cz, periodic); // build it once
  }

  for (int n=1; n <= nsteps; n++)        // repeat the loop nsteps times
  {
    if (!periodic) set_constant_boundary(Tnew, 1, Xdim, 1, Ydim, 1, Zdim, boundary);

    if (m == FTCS)
    {
      ftcs(Tnew, T, 1, Xdim, 1, Ydim, 1, Zdim, Cx, Cy, Cz, periodic);
    } else if (m == BE)
    {
      copy_dmatrix(A, constA, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim); // copy it every time :(
      becs(Tnew, T, 1, Xdim, 1, Ydim, 1, Zdim, A, periodic);
    }
    else if (m == CN)
    {
      copy_dmatrix(A, constA, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim); // copy it every time :(
      cn(Tnew, T, 1, Xdim, 1, Ydim, 1, Zdim, A, Cx, Cy, Cz, periodic);
    }

    copy_d3tensor(T, Tnew, 1, Xdim, 1, Ydim, 1, Zdim); // update T to the new values

    if (n%sample == 0)
    {
      printf("%s %d\n", methods[m], n);
      show_d3tensor("T", T, 1, Xdim, 1, Ydim, 1, Zdim);
      usleep(pause);
    }
  }

  free_dvector(x, 1, Xdim);
  free_dvector(y, 1, Ydim);
  free_dvector(z, 1, Zdim);
  if (m == BE || m == CN)
  {
    free_dmatrix(A, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
    free_dmatrix(constA, 1, Xdim*Ydim*Zdim, 1, Xdim*Ydim*Zdim);
  }
  free_d3tensor(T, 1, Xdim, 1, Ydim, 1, Zdim);
  free_d3tensor(Tnew, 1, Xdim, 1, Ydim, 1, Zdim);
  screen("[?25h"); // show cursor$
}

void populate_becs_matrix(double **A, long Xdim, long Ydim, long Zdim, 
                          double Cx, double Cy, double Cz, bool periodic)
{
  // Prepare matrix A from the the Backward Euler discretization using the Cx, Cy, Cz values

  // malloc doesn't guarentee zeroed memory so we zero it out first
  for (int i = 1; i <= Xdim*Ydim*Zdim; i++)
    for (int j = 1; j <= Xdim*Ydim*Zdim; j++)
      A[i][j] = 0;
  if (periodic)
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

void cn(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A, double Cx, double Cy, double Cz, bool periodic)
{
  if (periodic)
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
          long left, right;
          *B = src[i][j][k];
          left = (i == nrl) ? nrh : i-1; right = (i == nrh) ? nrl : i+1;
          *B += Cx*(src[left][j][k] + src[right][j][k] - 2*src[i][j][k]);
          left = (j == ncl) ? nch : j-1; right = (j == nch) ? ncl : j+1;
          *B += Cy*(src[i][left][k] + src[i][right][k] - 2*src[i][j][k]);
          left = (k == ndl) ? ndh : k-1; right = (k == ndh) ? ndl : k+1;
          *B += Cz*(src[i][j][left] + src[i][j][right] - 2*src[i][j][k]);
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
  } else {
    long X = nrh-nrl+1; // include constant boundaries, so they can influence other cells 
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
          *B = src[i][j][k] - 2*Cx*src[i][j][k] - 2*Cy*src[i][j][k] - 2*Cz*src[i][j][k];
          *B += (i==nrl) ? 0 : Cx*src[i-1][j][k];
          *B += (i==nrh) ? 0 : Cx*src[i+1][j][k];

          *B += (j==ncl) ? 0 : Cy*src[i][j-1][k];
          *B += (j==nch) ? 0 : Cy*src[i][j+1][k];

          *B += (k==ndl) ? 0 : Cz*src[i][j][k-1];
          *B += (k==ndh) ? 0 : Cz*src[i][j][k+1];
        }
    // show_dvector("b", b, 1, X*Y*Z);
  
    // solve Ax=b for x, using elimination
    gaussian_elimination(A, x, b, X*Y*Z);
    // show_dvector("x", x, 1, X*Y*Z);
  
    // unflatten x into dst 
    for (int i = nrl+1; i <= nrh-1; i++)
      for (int j = ncl+1; j <= nch-1; j++)
        for (int k = ndl+1; k <= ndh-1; k++)
          dst[i][j][k] = x[1+((i-nrl-1)*Y+(j-ncl-1))*Z+(k-ndl-1)];
    free_dvector(x, 1, X*Y*Z);
    free_dvector(b, 1, X*Y*Z);
  }
}

void becs(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double **A, bool periodic)
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

void ftcs(double ***dst, double ***src,
          long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
          double Cx, double Cy, double Cz, bool periodic)
{
 if (periodic)
  for (int i = nrl; i <= nrh; i++)            // calculate new values for all cells
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
      {
        long left, right;
        dst[i][j][k] = src[i][j][k];
        left = (i == nrl) ? nrh : i-1; right = (i == nrh) ? nrl : i+1;
        dst[i][j][k] += Cx*(src[left][j][k] + src[right][j][k] - 2*src[i][j][k]);
        left = (j == ncl) ? nch : j-1; right = (j == nch) ? ncl : j+1;
        dst[i][j][k] += Cy*(src[i][left][k] + src[i][right][k] - 2*src[i][j][k]);
        left = (k == ndl) ? ndh : k-1; right = (k == ndh) ? ndl : k+1;
        dst[i][j][k] += Cz*(src[i][j][left] + src[i][j][right] - 2*src[i][j][k]);
      }
  else
    for (int i = nrl+1; i <= nrh-1; i++)      // calculate new values for the non-boundaries
      for (int j = ncl+1; j <= nch-1; j++)
        for (int k = ndl+1; k <= ndh-1; k++)
          dst[i][j][k] = src[i][j][k] + Cx*(src[i-1][j][k] + src[i+1][j][k] - 2*src[i][j][k])
                                      + Cy*(src[i][j-1][k] + src[i][j+1][k] - 2*src[i][j][k])
                                      + Cz*(src[i][j][k-1] + src[i][j][k+1] - 2*src[i][j][k]);
}

void set_constant_boundary(double ***T,
                           long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
                           double value)
{
  // first row
  for (int j = ncl; j <= nch; j++)
    for (int k = ndl; k <= ndh; k++)
      T[nrl][j][k] = value;
  // middle rows
  for (int i = nrl+1; i < nrh; i++)
  {
    // first column, then middle columns, than last column
    for (int k = ndl; k <= ndh; k++) T[i][ncl][k] = value;
    for (int j = ncl+1; j < nch; j++) T[i][j][ndl] = T[i][j][ndh] = value;
    for (int k = ndl; k <= ndh; k++) T[i][nch][k] = value;
  }
  // last row
  for (int j = ncl; j <= nch; j++)
    for (int k = ndl; k <= ndh; k++)
      T[nrh][j][k] = value;
}

void set_initial_with_noise(double ***T,
                            long nrl, long nrh, long ncl, long nch, long ndl, long ndh,
                            double *x, double *y, double *z,
                            double(*init)(double,double,double), double max)
{
  // if the initial number is I, the result will be between (1-max)*I and (1+max)*I
  double A = 1 - max;
  double B = 2*max/RAND_MAX;
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
        T[i][j][k] = init(x[i], y[j], z[k])*(A + B*rand());
}

double gauss3(double x, double y, double z)
{
  return exp(-pow((5*x-2.5), 2))*exp(-pow((5*y-2.5), 2))*exp(-pow((5*z-2.5), 2));
}
