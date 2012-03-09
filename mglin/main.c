#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
int main(int argc, char **argv){
  FILE *outfile;
  double **f;
  int n=33;
  int ncycle=2;
  f = dmatrix(1,n,1,n);
  f[n/2+1][n/2+1]=-1000.0;
  H("mglin:");
  mglin(f,n,ncycle);
  H("fopen:");
  outfile = fopen("soln.dat", "w");
  fwrite(&f[1][1],sizeof(double),n*n,outfile);
  fclose(outfile);
  return 0;
}

void here(const char *s, const char *file, int line, const char *func){
  fprintf(stderr, "@ %s:%d in %s: %s\n", file, line, func, s);
}
