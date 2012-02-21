#ifndef UTILITIES_H
#define UTILITIES_H

#define DEBUG debug(__LINE__);
// #define DEBUG

void debug(int line);
void screen(const char *s);
void show_scalar(double x);
void show_dvector(const char*, double *v, long nl, long nh);
void show_dmatrix(const char* s, double **T, long nrl, long nrh, long ncl, long nch);
void copy_dmatrix(double **dst, double **src, long nrl, long nrh, long ncl, long nch);
void show_d3tensor(const char* s, double ***T, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void copy_d3tensor(double ***dst, long drl, long drh, long dcl, long dch, long ddl, long ddh,
                   double ***src, long srl, long srh, long scl, long sch, long sdl, long sdh);

// 3D array and bounds
typedef struct {
  double ***T;
  long nrl, nrh, ncl, nch, ndl, ndh;
} temp3D;

temp3D *create_temp3D(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void copy_temp3D(temp3D *d, const temp3D *s);
void show_temp3D(const char *s, const temp3D *t);
void free_temp3D(temp3D *t);

#endif
