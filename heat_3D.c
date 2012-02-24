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
  if (!p->quiet) screen("\033[2J\033[?25l"); // clear screen and hide cursor

  double *x, *y, *z;
  double **A, **constA;   // in case we use BE or CN
  int X = p->nx+2;
  int Y = p->ny+2;
  int Z = p->nz+2;
  x = dvector(1, X);
  y = dvector(1, Y);
  z = dvector(1, Z);
  for (int i = 1; i <= X; i++) x[i] = (i-1) / (double) (p->nx+1);
  for (int j = 1; j <= Y; j++) y[j] = (j-1) / (double) (p->ny+1);
  for (int k = 1; k <= Z; k++) z[k] = (k-1) / (double) (p->nz+1);

  t3D *t = create_t3D(1, X, 1, Y, 1, Z);
  t3D *tnew = create_t3D(1, X, 1, Y, 1, Z);
  
  if (!p->periodic) set_constant_boundary(p, t);
  set_initial_with_noise(p, t, x, y, z, init);

  if (!p->quiet) show_t3D("T", t);
  usleep(p->pause*2);

  if (p->method == BE || p->method == CN)
  {
    A = dmatrix(1, X*Y*Z, 1, X*Y*Z);
    constA = dmatrix(1, X*Y*Z, 1, X*Y*Z);
    if (p->method == CN) { Cx /= 2; Cy /= 2; Cz /= 2; }
    populate_becs_matrix(p, constA, X, Y, Z, Cx, Cy, Cz); // build it once
  }

  if (p->op) timer(true); // start timer outside loop
  for (int n=1; n <= p->nsteps; n++)        // repeat the loop nsteps times
  {
    if (p->method == FTCS)
    {
      ftcs(p, tnew, t, Cx, Cy, Cz);
    } else if (p->method == BE)
    {
      copy_dmatrix(A, constA, 1, X*Y*Z, 1, X*Y*Z); // copy it every time :(
      becs(p, tnew, t, A);
    }
    else if (p->method == CN)
    {
      copy_dmatrix(A, constA, 1, X*Y*Z, 1, X*Y*Z); // copy it every time :(
      cn(p, tnew, t, A, Cx, Cy, Cz);
    }

    copy_t3D(t, tnew);

    if (n%p->sample == 0)
    {
      if (p->op) { fprintf(p->op, "%d %.16e\n", n, timer(false)); }
      if (p->os) output_t3D(p->os, n, t);
      if (!p->quiet)
      {
        printf("%s %d\n", methodnames[p->method], n);
        show_t3D("T", t);
      }
      usleep(p->pause);
    }
  }

  free_dvector(x, 1, X);
  free_dvector(y, 1, Y);
  free_dvector(z, 1, Z);
  if (p->method == BE || p->method == CN)
  {
    free_dmatrix(A, 1, X*Y*Z, 1, X*Y*Z);
    free_dmatrix(constA, 1, X*Y*Z, 1, X*Y*Z);
  }
  free_t3D(t);
  free_t3D(tnew);
  if (!p->quiet) screen("\033[?25h"); // show cursor$
}

void populate_becs_matrix(const prefs3D *p, double **A, long X, long Y, long Z, 
                          double Cx, double Cy, double Cz)
{
  // Prepare matrix A from the the Backward Euler discretization using the Cx, Cy, Cz values

  // malloc doesn't guarentee zeroed memory so we zero it out first
  for (int i = 1; i <= X*Y*Z; i++)
    for (int j = 1; j <= X*Y*Z; j++)
      A[i][j] = 0;
  if (p->periodic)
  {
    for (int m = 1; m <= X*Y*Z; m++)
    {
      // Fill main diagonal
      A[m][m] = 2*(Cx+Cy+Cz)+1;
      int left, right;
  
      // Fill off diagonals
      left  = ((m-1)%Z == 0)   ? m+Z-1 : m-1;
      right = ((m-1)%Z == Z-1) ? m-Z+1 : m+1;
      A[m][left] = A[m][right] = -Cz;
  
      // Fill near bands
      left  = (((m-1)/Z)%Y == 0)   ? m+Z*(Y-1) : m-Z;
      right = (((m-1)/Z)%Y == Y-1) ? m-Z*(Y-1) : m+Z;
      A[m][left] = A[m][right] = -Cy;
  
      // Fill far bands
      left  = (((m-1)/(Z*Y)%X) == 0)   ? m+Z*Y*(X-1) : m-Z*Y;
      right = (((m-1)/(Z*Y)%X) == X-1) ? m-Z*Y*(X-1) : m+Z*Y;
      A[m][left] = A[m][right] = -Cx;
    }
  } else { // constant boundary
    for (int m = 1; m <= X*Y*Z; m++)
    {
      // if m is on any boundary, set diagonal to 1 and continue
      if ((m-1)%Z == 0 || (m-1)%Z == Z-1 ||
          ((m-1)/Z)%Y == 0 || ((m-1)/Z)%Y == Y-1 ||
          ((m-1)/(Z*Y)%X) == 0 || ((m-1)/(Z*Y)%X) == X-1)
      {
        A[m][m] = 1;
        continue;
      }
      // Fill main diagonal
      A[m][m] = 2*(Cx+Cy+Cz)+1;
      // Fill off diagonals
      A[m][m-1] = A[m][m+1] = -Cz;
      // Fill near bands
      A[m][m-Z] = A[m][m+Z] = -Cy;
      // Fill far bands
      A[m][m-Z*Y] = A[m][m+Z*Y] = -Cx;
    }
  }
  // show_dmatrix("A", A, 1, X*Y*Z, 1, X*Y*Z);
}

void cn(const prefs3D *p, t3D *d, t3D *s,
          double **A, double Cx, double Cy, double Cz)
{
  long X = s->nrh-s->nrl+1; // include boundaries
  long Y = s->nch-s->ncl+1;
  long Z = s->ndh-s->ndl+1;
  t3D *x = create_t3D(1, X, 1, Y, 1, Z);
  t3D *b = create_t3D(1, X, 1, Y, 1, Z);
  // flatten src into b
  for (int i = s->nrl; i <= s->nrh; i++)
    for (int j = s->ncl; j <= s->nch; j++)
      for (int k = s->ndl; k <= s->ndh; k++)
      {
        double *B = &b->T[i][j][k];
        *B = s->T[i][j][k];
        if (p->periodic) {
          long left, right;
          left = (i == s->nrl) ? s->nrh : i-1; right = (i == s->nrh) ? s->nrl : i+1;
          *B += Cx*(s->T[left][j][k] + s->T[right][j][k] - 2*s->T[i][j][k]);
          left = (j == s->ncl) ? s->nch : j-1; right = (j == s->nch) ? s->ncl : j+1;
          *B += Cy*(s->T[i][left][k] + s->T[i][right][k] - 2*s->T[i][j][k]);
          left = (k == s->ndl) ? s->ndh : k-1; right = (k == s->ndh) ? s->ndl : k+1;
          *B += Cz*(s->T[i][j][left] + s->T[i][j][right] - 2*s->T[i][j][k]);
        } else if (!(i==s->nrl || i==s->nrh || j==s->ncl ||j==s->nch || k==s->ndl || k==s->ndh)) {
          *B += - 2*Cx*s->T[i][j][k] - 2*Cy*s->T[i][j][k] - 2*Cz*s->T[i][j][k];
          *B += (i==s->nrl) ? 0 : Cx*s->T[i-1][j][k];
          *B += (i==s->nrh) ? 0 : Cx*s->T[i+1][j][k];
          *B += (j==s->ncl) ? 0 : Cy*s->T[i][j-1][k];
          *B += (j==s->nch) ? 0 : Cy*s->T[i][j+1][k];
          *B += (k==s->ndl) ? 0 : Cz*s->T[i][j][k-1];
          *B += (k==s->ndh) ? 0 : Cz*s->T[i][j][k+1];
        }
      }
  // solve Ax=b for x, using elimination (which expects 1-based 1D vectors)
  gaussian_elimination(A, (double *) x->T[1][1], (double *) b->T[1][1], X*Y*Z);
  
  // unflatten x into dst 
  copy_t3D(d, x);

  free_t3D(x);
  free_t3D(b);
}

void becs(const prefs3D *p, t3D *d, t3D *s,
          double **A)
{
  long X = s->nrh-s->nrl+1; // include boundaries
  long Y = s->nch-s->ncl+1;
  long Z = s->ndh-s->ndl+1;
  t3D *x = create_t3D(1, X, 1, Y, 1, Z);
  t3D *b = create_t3D(1, X, 1, Y, 1, Z);
  // flatten src into b
  copy_t3D(b, s);

  // solve Ax=b for x, using elimination (which expects 1-based 1D vectors)
  gaussian_elimination(A, (double *) x->T[1][1], (double *) b->T[1][1], X*Y*Z);

  // unflatten x into dst 
  copy_t3D(d, x);

  free_t3D(x);
  free_t3D(b);
}

void ftcs(const prefs3D *p, t3D *d, t3D *s,
          double Cx, double Cy, double Cz)
{
  for (int i = s->nrl; i <= s->nrh; i++)            // calculate new values for all cells
    for (int j = s->ncl; j <= s->nch; j++)
      for (int k = s->ndl; k <= s->ndh; k++)
      {
        d->T[i][j][k] = s->T[i][j][k];
        if (p->periodic) {
          long left, right;
          left = (i == s->nrl) ? s->nrh : i-1; right = (i == s->nrh) ? s->nrl : i+1;
          d->T[i][j][k] += Cx*(s->T[left][j][k] + s->T[right][j][k] - 2*s->T[i][j][k]);
          left = (j == s->ncl) ? s->nch : j-1; right = (j == s->nch) ? s->ncl : j+1;
          d->T[i][j][k] += Cy*(s->T[i][left][k] + s->T[i][right][k] - 2*s->T[i][j][k]);
          left = (k == s->ndl) ? s->ndh : k-1; right = (k == s->ndh) ? s->ndl : k+1;
          d->T[i][j][k] += Cz*(s->T[i][j][left] + s->T[i][j][right] - 2*s->T[i][j][k]);
        } else if (!(i==s->nrl || i==s->nrh || j==s->ncl ||j==s->nch || k==s->ndl || k==s->ndh)) {
            d->T[i][j][k] += Cx*(s->T[i-1][j][k] + s->T[i+1][j][k] - 2*s->T[i][j][k]) +
                            Cy*(s->T[i][j-1][k] + s->T[i][j+1][k] - 2*s->T[i][j][k]) +
                            Cz*(s->T[i][j][k-1] + s->T[i][j][k+1] - 2*s->T[i][j][k]);
        }
      }
}

void set_constant_boundary(const prefs3D *p, t3D *t)
{
  // first row
  for (int j = t->ncl; j <= t->nch; j++)
    for (int k = t->ndl; k <= t->ndh; k++)
      t->T[t->nrl][j][k] = p->boundary;
  // middle rows
  for (int i = t->nrl+1; i <= t->nrh-1; i++)
  {
    // first column, then middle columns, than last column
    for (int k = t->ndl; k <= t->ndh; k++) t->T[i][t->ncl][k] = p->boundary;
    for (int j = t->ncl+1; j <= t->nch-1; j++) t->T[i][j][t->ndl] = t->T[i][j][t->ndh] = p->boundary;
    for (int k = t->ndl; k <= t->ndh; k++) t->T[i][t->nch][k] = p->boundary;
  }
  // last row
  for (int j = t->ncl; j <= t->nch; j++)
    for (int k = t->ndl; k <= t->ndh; k++)
      t->T[t->nrh][j][k] = p->boundary;
}

void set_initial_with_noise(const prefs3D *p, t3D *t,
                            double *x, double *y, double *z,
                            double(*init)(double,double,double))
{
  // if the initial number is I, the result will be between (1-noise)*I and (1+noise)*I
  double A = 1 - p->noise;
  double B = 2*p->noise/RAND_MAX;
  long rl = t->nrl, rh = t->nrh, cl = t->ncl, ch = t->nch, dl = t->ndl, dh = t->ndh;
  if (!p->periodic) { rl++; rh--; cl++; ch--; dl++; dh--; }
  for (int i = rl; i <= rh; i++)
    for (int j = cl; j <= ch; j++)
      for (int k = dl; k <= dh; k++)
        t->T[i][j][k] = init(x[i], y[j], z[k])*(A + B*rand());
}

double gauss3(double x, double y, double z)
{
  return exp(-pow((5*x-2.5), 2))*exp(-pow((5*y-2.5), 2))*exp(-pow((5*z-2.5), 2));
}
