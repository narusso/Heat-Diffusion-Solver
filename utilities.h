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
void copy_d3tensor(double ***dst, double ***src, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

#endif
