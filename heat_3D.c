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
#include <assert.h> // assert

const char *methodnames[] = {
 [FTCS] = "Forward-Time Center-Space (explicit)",
   [BE] = "Backward Euler using Gaussian Elimination",
   [CN] = "Crank-Nicolson using Gaussian Elimination",
  [BEj] = "Backward Euler using Jacobi",
 [BEgs] = "Backward Euler using Gauss-Seidel",
[BEsor] = "Backward Euler using Successive over-relaxation",
};

void reset(int s)
{
  printf("\033[m");
  printf("\033[?25h");
  if (s != SIGSTOP)
    exit(EXIT_FAILURE);
}

d3D *create_d3D(long X, long Y, long Z)
{ // only uses 1-based vectors 
  d3D *disc;
  disc = (d3D *) malloc((size_t) sizeof(d3D));
  disc->xx = dvector(1, X);
  disc->yy = dvector(1, Y);
  disc->zz = dvector(1, Z);
  for (long i = 1; i <= X; i++) disc->xx[i] = (i-1) / (double) (X-1);
  for (long j = 1; j <= Y; j++) disc->yy[j] = (j-1) / (double) (Y-1);
  for (long k = 1; k <= Z; k++) disc->zz[k] = (k-1) / (double) (Z-1);
  return disc;
}

void free_d3D(d3D *disc)
{
  free_dvector(disc->xx, 1, disc->X);
  free_dvector(disc->yy, 1, disc->Y);
  free_dvector(disc->zz, 1, disc->Z);
  free(disc);
}

void solve(const prefs3D *p,
           double Cx, double Cy, double Cz,
           double(*init)(double,double,double))
{
  signal(SIGFPE, reset);  // works
  signal(SIGINT, reset);  // works
  signal(SIGPIPE, reset); // doesn't reset cursor
  if (!p->quiet) screen("\033[2J\033[?25l"); // clear screen and hide cursor

  double **A, **constA;   // in case we use BE or CN
  long X = p->nx+2;
  long Y = p->ny+2;
  long Z = p->nz+2;
  d3D *disc = create_d3D(X, Y, Z);

  t3D *t = create_t3D(1, X, 1, Y, 1, Z);
  t3D *tnew = create_t3D(1, X, 1, Y, 1, Z);
  
  if (!p->periodic) set_constant_boundary(p, t);
  set_initial_with_noise(p, t, disc, init);

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
      ftcs(p, tnew, t, Cx, Cy, Cz, disc);
    } else if (p->method == BE)
    {
      copy_dmatrix(A, constA, 1, X*Y*Z, 1, X*Y*Z); // copy it every time :(
      becs(p, tnew, t, A, disc);
    }
    else if (p->method == CN)
    {
      copy_dmatrix(A, constA, 1, X*Y*Z, 1, X*Y*Z); // copy it every time :(
      cn(p, tnew, t, A, Cx, Cy, Cz, disc);
    }
    else if (p->method == BEj)
    {
      bej(p, tnew, t, Cx, Cy, Cz, disc);
    }
    else if (p->method == BEgs)
    {
      begs(p, tnew, t, Cx, Cy, Cz, disc);
    }
    else if (p->method == BEsor)
    {
      besor(p, tnew, t, Cx, Cy, Cz, disc);
    }

    copy_t3D(t, tnew);

    if (n%p->sample == 0)
    {
      if (p->op) { fprintf(p->op, "%d %.3e\n", n, timer(false)); }
      if (p->os) output_t3D(p->os, n, t);
      if (!p->quiet)
      {
        printf("%s %d\n", methodnames[p->method], n);
        show_t3D("T", t);
      }
      usleep(p->pause);
    }
  }

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
  for (long i = 1; i <= X*Y*Z; i++)
    for (long j = 1; j <= X*Y*Z; j++)
      A[i][j] = 0;
  if (p->periodic)
  {
    for (long m = 1; m <= X*Y*Z; m++)
    {
      // Fill main diagonal
      A[m][m] = 2*(Cx+Cy+Cz)+1;
      long left, right;
  
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
    for (long m = 1; m <= X*Y*Z; m++)
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

void bej(const prefs3D *p, t3D *d, t3D *s,
         double Cx, double Cy, double Cz, const d3D *disc)
{
  // determine values for d->T from s->T according to prefs from p
  // and constants Cx, Cy, Cz
  // pretend all sorts of things for now (eg. boundary = 0)
  // pretend d, s, and temp t3D all use same indexing
  // pretend temp->T is zeroed out
  copy_t3D(d, s); // set first guess to last solution
  t3D *temp = create_t3D(d->nrl, d->nrh, d->ncl, d->nch, d->ndl, d->ndh);
  double ***x = d->T;
  double ***xnew = temp->T;
  double ***xold = s->T;
  int MAX_ITER = 2000;
  long elements = (s->nrh-s->nrl+1)*(s->nch-s->ncl+1)*(s->ndh-s->ndl+1);
  for (long m = 0; m < MAX_ITER; m++)
  {
    double diff = 0;
    for (long i = s->nrl+1; i <= s->nrh-1; i++)
      for (long j = s->ncl+1; j <= s->nch-1; j++)
        for (long k = s->ndl+1; k <= s->ndh-1; k++)
        {
          xnew[i][j][k] = Cx/(2*Cx+1)*(x[i-1][j][k] + x[i+1][j][k])
                        + Cy/(2*Cy+1)*(x[i][j-1][k] + x[i][j+1][k])
                        + Cz/(2*Cz+1)*(x[i][j][k-1] + x[i][j][k+1])
                        + 1/(2*Cx+2*Cy+2*Cz+1) * xold[i][j][k];
          diff += fabs(xnew[i][j][k] - x[i][j][k]);
        }
    if (diff/elements < 1.e-15) break;
    double ***t = x; x = xnew; xnew = t;
  }
  if (p->source)
    for (long i = d->nrl+1; i <= d->nrh-1; i++)
      for (long j = d->ncl+1; j <= d->nch-1; j++)
        for (long k = d->ndl+1; k <= d->ndh-1; k++)
          d->T[i][j][k] += p->source(i, j, k, disc);
  free_t3D(temp);
}

void begs(const prefs3D *p, t3D *d, t3D *s,
          double Cx, double Cy, double Cz, const d3D *disc)
{
  // determine values for d->T from s->T according to prefs from p
  // and constants Cx, Cy, Cz
  // pretend all sorts of things for now (eg. boundary = 0)
  // pretend d, s, and temp t3D all use same indexing
  // pretend temp->T is zeroed out
  copy_t3D(d, s); // set first guess to last solution
  double ***x = d->T;
  double ***xold = s->T;
  int MAX_ITER = 2000;
  long elements = (s->nrh-s->nrl+1)*(s->nch-s->ncl+1)*(s->ndh-s->ndl+1);
  for (long m = 0; m < MAX_ITER; m++)
  {
    double diff = 0;
    for (long i = s->nrl+1; i <= s->nrh-1; i++)
      for (long j = s->ncl+1; j <= s->nch-1; j++)
        for (long k = s->ndl+1; k <= s->ndh-1; k++)
        {
          double t = Cx/(2*Cx+1)*(x[i-1][j][k] + x[i+1][j][k])
                   + Cy/(2*Cy+1)*(x[i][j-1][k] + x[i][j+1][k])
                   + Cz/(2*Cz+1)*(x[i][j][k-1] + x[i][j][k+1])
                   + 1/(2*Cx+2*Cy+2*Cz+1) * xold[i][j][k];
          diff += fabs(t - x[i][j][k]);
          x[i][j][k] = t;
        }
    if (diff/elements < 1.e-15) break;
  }
  if (p->source)
    for (long i = d->nrl+1; i <= d->nrh-1; i++)
      for (long j = d->ncl+1; j <= d->nch-1; j++)
        for (long k = d->ndl+1; k <= d->ndh-1; k++)
          d->T[i][j][k] += p->source(i, j, k, disc);
}

void besor(const prefs3D *p, t3D *d, t3D *s,
           double Cx, double Cy, double Cz, const d3D *disc)
{
  double w = 1.65; // move into prefs
  // determine values for d->T from s->T according to prefs from p
  // and constants Cx, Cy, Cz
  // pretend all sorts of things for now (eg. boundary = 0)
  // pretend d, s, and temp t3D all use same indexing
  // pretend temp->T is zeroed out
  copy_t3D(d, s); // set first guess to last solution
  double ***x = d->T;
  double ***xold = s->T;
  int MAX_ITER = 2000;
  long elements = (s->nrh-s->nrl+1)*(s->nch-s->ncl+1)*(s->ndh-s->ndl+1);
  for (long m = 0; m < MAX_ITER; m++)
  {
    double diff = 0;
    for (long i = s->nrl+1; i <= s->nrh-1; i++)
      for (long j = s->ncl+1; j <= s->nch-1; j++)
        for (long k = s->ndl+1; k <= s->ndh-1; k++)
        {
          double t = (1-w)*x[i][j][k] + w*Cx/(2*Cx+1)*(x[i-1][j][k] + x[i+1][j][k])
                   + w*Cy/(2*Cy+1)*(x[i][j-1][k] + x[i][j+1][k])
                   + w*Cz/(2*Cz+1)*(x[i][j][k-1] + x[i][j][k+1])
                   + w/(2*Cx+2*Cy+2*Cz+1) * xold[i][j][k];
          diff += fabs(t - x[i][j][k]);
          x[i][j][k] = t;
        }
    if (diff/elements < 1.e-15) break;
  }
  if (p->source)
    for (long i = d->nrl+1; i <= d->nrh-1; i++)
      for (long j = d->ncl+1; j <= d->nch-1; j++)
        for (long k = d->ndl+1; k <= d->ndh-1; k++)
          d->T[i][j][k] += p->source(i, j, k, disc);
}

void cn(const prefs3D *p, t3D *d, t3D *s,
          double **A, double Cx, double Cy, double Cz, const d3D *disc)
{
  long X = s->nrh-s->nrl+1; // include boundaries
  long Y = s->nch-s->ncl+1;
  long Z = s->ndh-s->ndl+1;
  t3D *x = create_t3D(1, X, 1, Y, 1, Z);
  t3D *b = create_t3D(1, X, 1, Y, 1, Z);
  // flatten src into b
  for (long i = s->nrl; i <= s->nrh; i++)
    for (long j = s->ncl; j <= s->nch; j++)
      for (long k = s->ndl; k <= s->ndh; k++)
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

  if (p->source)
    for (long i = d->nrl+1; i <= d->nrh-1; i++)
      for (long j = d->ncl+1; j <= d->nch-1; j++)
        for (long k = d->ndl+1; k <= d->ndh-1; k++)
          d->T[i][j][k] += p->source(i, j, k, disc);

  free_t3D(x);
  free_t3D(b);
}

void becs(const prefs3D *p, t3D *d, t3D *s,
          double **A, const d3D *disc)
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

  if (p->source)
    for (long i = d->nrl+1; i <= d->nrh-1; i++)
      for (long j = d->ncl+1; j <= d->nch-1; j++)
        for (long k = d->ndl+1; k <= d->ndh-1; k++)
          d->T[i][j][k] += p->source(i, j, k, disc);

  free_t3D(x);
  free_t3D(b);
}

void ftcs(const prefs3D *p, t3D *d, t3D *s,
          double Cx, double Cy, double Cz, const d3D *disc)
{
  for (long i = s->nrl; i <= s->nrh; i++)            // calculate new values for all cells
    for (long j = s->ncl; j <= s->nch; j++)
      for (long k = s->ndl; k <= s->ndh; k++)
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
          if (p->source) d->T[i][j][k] += p->source(i, j, k, disc);
        } else if (!(i==s->nrl || i==s->nrh || j==s->ncl ||j==s->nch || k==s->ndl || k==s->ndh)) {
          d->T[i][j][k] += Cx*(s->T[i-1][j][k] + s->T[i+1][j][k] - 2*s->T[i][j][k]) +
                           Cy*(s->T[i][j-1][k] + s->T[i][j+1][k] - 2*s->T[i][j][k]) +
                           Cz*(s->T[i][j][k-1] + s->T[i][j][k+1] - 2*s->T[i][j][k]);
          if (p->source) d->T[i][j][k] += p->source(i, j, k, disc);
        }
      }
}

void set_constant_boundary(const prefs3D *p, t3D *t)
{
  // first row
  for (long j = t->ncl; j <= t->nch; j++)
    for (long k = t->ndl; k <= t->ndh; k++)
      t->T[t->nrl][j][k] = p->boundary;
  // middle rows
  for (long i = t->nrl+1; i <= t->nrh-1; i++)
  {
    // first column, then middle columns, than last column
    for (long k = t->ndl; k <= t->ndh; k++) t->T[i][t->ncl][k] = p->boundary;
    for (long j = t->ncl+1; j <= t->nch-1; j++) t->T[i][j][t->ndl] = t->T[i][j][t->ndh] = p->boundary;
    for (long k = t->ndl; k <= t->ndh; k++) t->T[i][t->nch][k] = p->boundary;
  }
  // last row
  for (long j = t->ncl; j <= t->nch; j++)
    for (long k = t->ndl; k <= t->ndh; k++)
      t->T[t->nrh][j][k] = p->boundary;
}

void set_initial_with_noise(const prefs3D *p, t3D *t,
                            const d3D *disc,
                            double(*init)(double,double,double))
{
  // if the initial number is I, the result will be between (1-noise)*I and (1+noise)*I
  double A = 1 - p->noise;
  double B = 2*p->noise/RAND_MAX;
  long rl = t->nrl, rh = t->nrh, cl = t->ncl, ch = t->nch, dl = t->ndl, dh = t->ndh;
  if (!p->periodic) { rl++; rh--; cl++; ch--; dl++; dh--; }
  for (long i = rl; i <= rh; i++)
    for (long j = cl; j <= ch; j++)
      for (long k = dl; k <= dh; k++)
        t->T[i][j][k] = init(disc->xx[i], disc->yy[j], disc->zz[k])*(A + B*rand());
}

double gauss3(double x, double y, double z)
{
  assert(x <= 1); assert(x >= 0); assert(y <= 1); assert(y >= 0); assert(z <= 1); assert(z >= 0);
  return exp(-pow((5*x-2.5), 2))*exp(-pow((5*y-2.5), 2))*exp(-pow((5*z-2.5), 2));
}

double plane_source(long i, long j, long k, const d3D *disc)
{
  double x = disc->xx[i];
  double y = disc->yy[j];
  double z = disc->zz[k];
  assert(x <= 1); assert(x >= 0); assert(y <= 1); assert(y >= 0); assert(z <= 1); assert(z >= 0);
  if (fabs(x-y) < 1e-4) return .00000001;
  return 0;
}
