#include "utilities.h"
#include <stdio.h> // printf

void debug(int line) { printf("%d\n", line); }
void screen(const char *s) { printf(s); }
// void screen(const char *s) {}

void show_dvector(const char* s, double *v, long nl, long nh)
{
  printf("%s:\n", s);
  for (long i = nl; i <= nh; i++)
  {
    show_scalar(v[i]);
    printf(i%8 ? " " : "\n");
  }
  printf("\n");
}

void show_dmatrix(const char* s, double **T, long nrl, long nrh, long ncl, long nch)
{
  screen("\033[1H"); // move to top left corner of the screen
  printf("%s:\n", s);
  for (long i = nrl; i <= nrh; i++)
  {
    printf("%3ld: ", i);
    for (long j = ncl; j <= nch; j++)
      show_scalar(T[i][j]);
    printf("\n");
  }
  printf("\n");
}

void show_d3tensor(const char* s, double ***T, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  screen("\033[1H"); // move to top left corner of the screen
  printf("%s:\n", s);
  for (long k = ndl; k <= ndh; k++)
  {
    printf("k=%3ld\n", k);
    for (long i = nrl; i <= nrh; i++)
    {
      printf("i=%3ld: ", i);
      for (long j = ncl; j <= nch; j++)
        show_scalar(T[i][j][k]);
      printf("\n");
    }
    printf("\n");
  }
}

void copy_d3tensor(double ***dst, double ***src, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      for (int k = ndl; k <= ndh; k++)
        dst[i][j][k] = src[i][j][k];
}

void copy_dmatrix(double **dst, double **src, long nrl, long nrh, long ncl, long nch)
{
  for (int i = nrl; i <= nrh; i++)
    for (int j = ncl; j <= nch; j++)
      dst[i][j] = src[i][j];
}

void show_scalar(double x)
{
  if (x == 0) printf("\033[33m");
  printf("%8.4f", x);
  if (x == 0) printf("\033[m");
}
