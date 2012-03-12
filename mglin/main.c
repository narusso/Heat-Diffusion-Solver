#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
int main(int argc, char **argv){
  FILE *outfile;
  double **f;
  int n=9;
  int ncycle=2;
  int nsteps=2;
  int show=1;
  int t; // current timestep
  f = dmatrix(1,n,1,n);
  int i, j;
  for (i=2;i<n;i++)
    for (j=2;j<n;j++)
      f[i][j]=30.0;
  f[n/2+1][n/2+1]=80.0;
  H("fopen:");
  outfile = fopen("soln.dat", "w");
  H("mglin:");
  for (t=0;t<nsteps;t++)
  {
    // Each time we call mglin, f is the solution at the latest timestep
    // for BE, we solve Ax=b for x, given
    // A, which depends on dt, alpha, dx, and
    // b, the temps at previous timestep
    // A will be hardcoded into slvsml and relax
    // source term will be there too
    if (t%show==0) fwrite(&f[1][1],sizeof(double),n*n,outfile);
    mglin(f,n,ncycle);
  }
  fclose(outfile);
  return 0;
}

void here(const char *s, const char *file, int line, const char *func){
  fprintf(stderr, "@ %s:%d in %s: %s\n", file, line, func, s);
}
