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

void mgsolve(const prefs3D *p)
{
  const int ncycle=3, NPRE=2, NPOST=2;
  signal(SIGFPE, reset);  // works
  signal(SIGINT, reset);  // works
  signal(SIGPIPE, reset); // doesn't reset cursor
  if (!p->quiet) screen("\033[2J\033[?25l"); // clear screen and hide cursor

  // Create temperature matrix of appropriate size
  long X = p->nx+2;
  long Y = p->ny+2;
  long Z = p->nz+2;
  t3D *t = create_t3D(1, X, 1, Y, 1, Z); // corresponds to u & n in mglin

  // Initialize matrix
  set_constant_boundary(p, t);
  d3D *disc = create_d3D(X, Y, Z);
  set_initial_with_noise(p, t, disc, p->init);

  // Show temps
  if (!p->quiet) show_t3D("T", t);
  usleep(p->pause*2);

  // Find number of grid levels needed
  int ng = 0;
  int nn = X; // we require nx==ny==nz
  while (nn >>= 1) ng++;
  if (X != 1+(1L << ng)) nrerror("nx+1 must be a power of 2 for multigrid.");

  // Create arrays of pointers to data at each grid level
  double ***solution[ng+1]; // +1 to allow for 1-based addressing?
  double ***rhs[ng+1];
  double ***res[ng+1];
  double ***destination[ng+1]; // not sure about rho vs rhs, so naming this generically for now

  int ngrid = ng-1; // ng is finest (original), ng-1 is 2nd finest, 1 is coarsest
  nn = X/2+1; // new total number of points

  // About here I need a master loop to handle multiple timesteps
  if (p->op) timer(true); // start timer outside loop
  for (int n=1; n <= p->nsteps; n++)        // repeat the loop nsteps times
  {
    // Restrict solution to next (coarser) grid
    destination[ngrid] = d3tensor(1, nn, 1, nn, 1, nn);
    rstrct(destination[ngrid], t->T, nn);
  
    while (nn > 3) {
      nn = nn/2+1;
      destination[--ngrid] = d3tensor(1, nn, 1, nn, 1, nn);
      rstrct(destination[ngrid], destination[ngrid+1], nn);
    }
  
    // Solve directly at coarsest level
    nn=3;
    solution[1] = d3tensor(1,nn,1,nn,1,nn);
    rhs[1]= d3tensor(1,nn,1,nn,1,nn); // why?
    slvsml(p, solution[1], destination[1]);
    free_d3tensor(destination[1],1,nn,1,nn,1,nn);
    ngrid=ng;
  
    for (int j=2;j<=ngrid;j++) {        // loop over coarse to fine, starting at level 2
      fprintf(stderr, "at grid level %d\n",j);
      nn=2*nn-1;
      solution[j]=d3tensor(1,nn,1,nn,1,nn);     // setup grids for lhs,rhs, and residual
      rhs[j]=d3tensor(1,nn,1,nn,1,nn);
      res[j]=d3tensor(1,nn,1,nn,1,nn);
      interp(solution[j],solution[j-1],nn);
      
      // destination contains rhs except on fine grid where it is in t->T
      copy(rhs[j],(j != ngrid ? destination[j] : t->T),nn);
  
      // v-cycle at current grid level
      for (int jcycle=1;jcycle<=ncycle;jcycle++) {
        fprintf(stderr, "vcycle number %d\n", jcycle);
      
        // nf is # points on finest grid for current v-sweep
        int nf=nn, jj;
        for (jj=j;jj>=2;jj--) {
          for (int jpre=1;jpre<=NPRE;jpre++)     // NPRE g-s sweeps on the finest (relatively) scale
            relax(p, solution[jj],rhs[jj],nf);
          resid(p, res[jj],solution[jj],rhs[jj],nf); // compute res on finest scale, store in res
          nf=nf/2+1;                             // next coarsest scale
          rstrct(rhs[jj-1],res[jj],nf);          // restrict residuals as rhs of next coarsest scale
          fill0(solution[jj-1],nf);              // set the initial solution guess to zero
        }
        slvsml(p, solution[1],rhs[1]);           // solve the small problem exactly
        nf=3;                                    // fine scale now n=3
        for (jj=2;jj<=j;jj++) {                  // work way back up to current finest grid
          nf=2*nf-1;                             // next finest scale
          addint(solution[jj],solution[jj-1],res[jj],nf);  // inter error and add to previous solution guess
          for (int jpost=1;jpost<=NPOST;jpost++) // do NPOST g-s sweeps
            relax(p, solution[jj],rhs[jj],nf);
        }
      }
    }

    copy(t->T, solution[ngrid], n);
    fprintf(stderr, "%20.9f\n", t->T[4][4][4]);
  
    if (n%p->sample == 0)
    {
      if (p->op) { fprintf(p->op, "%d %.3e\n", n, timer(false)); }
      if (p->os) output_t3D(p->os, n, t);
      if (!p->quiet)
      {
        printf("MG %d\n", n);
        show_t3D("T", t);
      }
      usleep(p->pause);
    }
  }

  // free resources
  free_t3D(t);
  free_d3D(disc);
  if (!p->quiet) screen("\033[?25h"); // show cursor

  return;
}

void solve(const prefs3D *p)
{
  signal(SIGFPE, reset);  // works
  signal(SIGINT, reset);  // works
  signal(SIGPIPE, reset); // doesn't reset cursor
  if (!p->quiet) screen("\033[2J\033[?25l"); // clear screen and hide cursor

  // derived constants
  double dx, dy, dz, Cx, Cy, Cz;
  dx = p->LX/(p->nx+1); dy = p->LY/(p->ny+1); dz = p->LZ/(p->nz+1);
  Cx = p->alpha*p->dt/(dx*dx);
  Cy = p->alpha*p->dt/(dy*dy);
  Cz = p->alpha*p->dt/(dz*dz);

  double **A, **constA;   // in case we use BE or CN
  long X = p->nx+2;
  long Y = p->ny+2;
  long Z = p->nz+2;
  d3D *disc = create_d3D(X, Y, Z);

  t3D *t = create_t3D(1, X, 1, Y, 1, Z);
  t3D *tnew = create_t3D(1, X, 1, Y, 1, Z);

  if (!p->periodic) set_constant_boundary(p, t);
  set_initial_with_noise(p, t, disc, p->init);

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
  free_d3D(disc);
  if (!p->quiet) screen("\033[?25h"); // show cursor
}

void populate_becs_matrix(const prefs3D *p, double **A, long X, long Y, long Z,
                          double Cx, double Cy, double Cz)
{
  // Prepare matrix A from the Backward Euler discretization using the Cx, Cy, Cz values

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
          xnew[i][j][k] = Cx/(2*Cx+2*Cy+2*Cz+1)*(x[i-1][j][k] + x[i+1][j][k])
                        + Cy/(2*Cx+2*Cy+2*Cz+1)*(x[i][j-1][k] + x[i][j+1][k])
                        + Cz/(2*Cx+2*Cy+2*Cz+1)*(x[i][j][k-1] + x[i][j][k+1])
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
          double t = Cx/(2*Cx+2*Cy+2*Cz+1)*(x[i-1][j][k] + x[i+1][j][k])
                   + Cy/(2*Cx+2*Cy+2*Cz+1)*(x[i][j-1][k] + x[i][j+1][k])
                   + Cz/(2*Cz+2*Cy+2*Cz+1)*(x[i][j][k-1] + x[i][j][k+1])
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
  // number of iterations is 0 when w=0, as then it just copies the old value.
  // number of iterations is low when w=1, which turns this into gauss-seidel.
  // shouldn't there be some other good value?
  double w = p->w;
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
          double t = (1-w)*x[i][j][k]
                   + w*Cx/(2*Cx+2*Cy+2*Cz+1)*(x[i-1][j][k] + x[i+1][j][k])
                   + w*Cy/(2*Cx+2*Cy+2*Cz+1)*(x[i][j-1][k] + x[i][j+1][k])
                   + w*Cz/(2*Cx+2*Cy+2*Cz+1)*(x[i][j][k-1] + x[i][j][k+1])
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
  if (fabs(x-y) < 1e-4) return .0001;
  return 0;
}

void rstrct(double ***dst, double ***src, long nh) // coarse nrl==ncl==ndl==1, nrh==nch==ndh==nh
{
  // for every point in coarse grid, if it's boundary, copy directly, else use half-weighting.
  // don't bother incrementing indices into fine grid, as they are easily calculated.
  long pf, qf, rf; // 3D fine indices: 1 to nh by 2's (eg. 1 3 5 7 9)
  long pc, qc, rc; // 3D coarse indices: 1 to nh*2-1  (eg. 1 2 3 4 5)
  for (pc=1; pc <= nh; pc++)
    for (qc=1; qc <= nh; qc++)
      for (rc=1; rc <= nh; rc++)
      {
        pf = pc*2-1; qf = qc*2-1; rf = rc*2-1;
        if (pc==1 || pc==nh || qc==1 || qc == nh || rc==1 || rc == nh)
          dst[pc][qc][rc] = src[pf][qf][rf];
        else
          dst[pc][qc][rc] = src[pf][qf][rf]/2.0 +
                           (src[pf][qf][rf-1] + src[pf][qf][rf+1] +
                            src[pf][qf-1][rf] + src[pf][qf+1][rf] +
                            src[pf-1][qf][rf] + src[pf+1][qf][rf])/12.0;
        /*printf("d[%ld][%ld][%ld]: %f s[%ld][%ld][%ld]: %f\n",
               pc,qc,rc,dst[pc][qc][rc],
               pf,qf,rf,src[pf][qf][rf]);*/
      }
}

void interp(double ***dst, double ***src, long nh) // fine nrl==ncl==ndl==1, nrh==nch==ndh==nh
{
  long pf, qf, rf; // 3D fine indices: 1 to nh by 2's (eg. 1 3 5 7 9)
  for (pf=1; pf <= nh; pf+=2) // odd, odd, odd
    for (qf=1; qf <= nh; qf+=2)
      for (rf=1; rf <= nh; rf+=2)
        dst[pf][qf][rf] = src[(pf+1)/2][(qf+1)/2][(rf+1)/2];
  for (pf=2; pf <= nh; pf+=2) // even, odd, odd
    for (qf=1; qf <= nh; qf+=2)
      for (rf=1; rf <= nh; rf+=2)
        dst[pf][qf][rf] = 0.5*(dst[pf-1][qf][rf]+dst[pf+1][qf][rf]);
  for (pf=1; pf <= nh; pf+=1) // all, even, odd
    for (qf=2; qf <= nh; qf+=2)
      for (rf=1; rf <= nh; rf+=2)
        dst[pf][qf][rf] = 0.5*(dst[pf][qf-1][rf]+dst[pf][qf+1][rf]);
  for (pf=1; pf <= nh; pf+=1) // all, all, even
    for (qf=1; qf <= nh; qf+=1)
      for (rf=2; rf <= nh; rf+=2)
        dst[pf][qf][rf] = 0.5*(dst[pf][qf][rf-1]+dst[pf][qf][rf+1]);
}

void slvsml(const prefs3D *p, double ***solution, double ***rhs)
{
  double dx=0.5;
  double C = p->alpha*p->dt/(dx*dx);
  fill0(solution,3);
  solution[2][2][2] = rhs[2][2][2]/(1+6*C);
}

void fill0(double ***m, long n)
{
  for (int p=1; p<=n; p++)
    for (int q=1; q<=n; q++)
      for (int r=1; r<=n; r++)
        m[p][q][r] = 0;
}

void copy(double ***dst, double ***src, long nh)
{
  for (long p=1; p<=nh; p++)
    for (long q=1; q<=nh; q++)
      for (long r=1; r<=nh; r++)
        dst[p][q][r]=src[p][q][r];
}

void stepby(int *i, int *j, int *k, int d, int n)
{
  *k += d;
  if (*k <= n) return;
  *k -= n;
  *j += 1;
  if (*j <= n) return;
  *j -= n;
  *i += 1;
}

void relax(const prefs3D *p, double ***sol, double ***rhs, int n)
{
  double dx = 1.0/(n-1);
  double C = p->alpha*p->dt/(dx*dx);
  for (int start=1; start <= 2; start++)
    for (int i=1, j=1, k=start; i<=n && j<=n && k<=n; stepby(&i, &j, &k, 2, n))
    {
      if (i==1 || i==n || j==1 || j==n || k==1 || k == n) continue; // skip boundaries
      sol[i][j][k] = C/(4*C+1)*(sol[i+1][j][k]+sol[i-1][j][k]+
                                sol[i][j+1][k]+sol[i][j-1][k]+
                                sol[i][j][k+1]+sol[i][j][k-1])
                   + 1/(4*C+1)*rhs[i][j][k];
    }
}

void resid(const prefs3D *p, double ***res, double ***sol, double ***rhs, int nf)
{
  double dx = 1.0/(nf-1);
  double C = p->alpha*p->dt/(dx*dx);
  for (int p=1; p <= nf; p++)
    for (int q=1; q <= nf; q++)
      for (int r=1; r <= nf; r++)
      {
        if (p==1 || p==nf || q==1 || q == nf || r==1 || r == nf)
          res[p][q][r] = 0.0;
        else
          // I've no reason to think this is correct:
          res[p][q][r] = (sol[p][q][r-1] + sol[p][q][r+1] +
                          sol[p][q-1][r] + sol[p][q+1][r] +
                          sol[p-1][q][r] + sol[p+1][q][r] +
                          sol[p][q][r]*(1+6*C))/(-dx*dx) +
                          rhs[p][q][r];
      }
}

void addint(double ***solcur, double ***solprev, double ***res, int nf)
{
  interp(solcur, solprev, nf);
  for (long p=1; p<=nf; p++)
    for (long q=1; q<=nf; q++)
      for (long r=1; r<=nf; r++)
        solcur[p][q][r]=res[p][q][r];
}
