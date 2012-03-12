#ifndef MG_H
#define MG_H

#define NPRE  2
#define NPOST 2
#define NGMAX 15

void mglin(double **u, int n, int ncycle);
void addint(double **uf, double **uc, double **res, int nf);
void copy(double **aout, double **ain, int n);
void fill0(double **u, int n);
void interp(double **uf, double **uc, int nf);
void relax(double **u, double **rhs, int n);
void resid(double **res, double **u, double **rhs, int n);
void rstrct(double **uc, double **uf, int nc);
void slvsml(double **u, double **rhs);
void here(const char *s, const char *file, int line, const char *func);
#define H(S) here(S,__FILE__,__LINE__,__func__);

#endif
