#include "utilities.h"
#include "nrutil.h"
#include <stdlib.h> // malloc and free
#include <stdio.h>  // printf
#include <assert.h> // assert
#include <time.h>   // CLOCKS_PER_SEC
#include <sys/times.h> // times
#include <unistd.h> // sysconf

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

void copy_d3tensor(double ***dst, long drl, long drh, long dcl, long dch, long ddl, long ddh,
                   double ***src, long srl, long srh, long scl, long sch, long sdl, long sdh)
{
  long oi = drl - srl;
  long oj = dcl - scl;
  long ok = ddl - sdl;

  for (int i = srl; i <= srh; i++)
    for (int j = scl; j <= sch; j++)
      for (int k = sdl; k <= sdh; k++)
        dst[i+oi][j+oj][k+ok] = src[i][j][k];
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

t3D *create_t3D(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  t3D *t;
  t = (t3D *) malloc((size_t) sizeof(t3D));
  if (!t) nrerror("allocation failure in create_t3D()");
  t->T = d3tensor(nrl, nrh, ncl, nch, ndl, ndh);
  t->nrl = nrl; t->nrh = nrh;
  t->ncl = ncl; t->nch = nch;
  t->ndl = ndl; t->ndh = ndh;
  return t;
}

void copy_t3D(t3D *d, const t3D *s)
{
  assert(d->nrh - d->nrl == s->nrh - s->nrl);
  assert(d->nch - d->ncl == s->nch - s->ncl);
  assert(d->ndh - d->ndl == s->ndh - s->ndl);
  copy_d3tensor(d->T, d->nrl, d->nrh, d->ncl, d->nch, d->ndl, d->ndh,
                s->T, s->nrl, s->nrh, s->ncl, s->nch, s->ndl, s->ndh);
}

void show_t3D(const char *s, const t3D *t)
{
  show_d3tensor(s, t->T, t->nrl, t->nrh, t->ncl, t->nch, t->ndl, t->ndh);
}

void output_t3D(FILE *stream, int n, const t3D *t)
{
  for (long i = t->nrl; i <= t->nrh; i++)
    for (long j = t->ncl; j <= t->nch; j++)
      for (long k = t->ndl; k <= t->ndh; k++)
        fprintf(stream, "%d %ld %ld %ld %.16e\n", n, i, j, k, t->T[i][j][k]);
}

void free_t3D(t3D *t)
{
  free_d3tensor(t->T, t->nrl, t->nrh, t->ncl, t->nch, t->ndl, t->ndh);
  t->T = NULL;
  free(t);
}

double timer(bool mode)
{
  // true: start the timer
  // false: return time in seconds since timer was last started
  static clock_t initial_time;
  static long ticks_per_second = 0;
  if (ticks_per_second==0) ticks_per_second = sysconf(_SC_CLK_TCK);
  struct tms buf;
  times(&buf);
  if (mode) initial_time = buf.tms_utime+buf.tms_stime;
  return (buf.tms_utime+buf.tms_stime-initial_time) / (double) ticks_per_second;
}
