#include "mg.h"
void slvsml(double **u, double **rhs)
/* 
   Solution of the model problem on the coarsest grid, where h = 1
   2 . The right-hand side is input
   in rhs[1..3][1..3] and the solution is returned in u[1..3][1..3].
*/
{
  H("");
  float alpha = 1.1234e-4; // diffusivity of copper in m^2/s
  float dt = 0.00001;      // timestep chosen arbitrarily
  float dx = 0.5;          // domain is [0,1]; 2 intervals at this coarseness
  float C = alpha*dt/(dx*dx);
  void fill0(double **u, int n);
  fill0(u,3);
  u[2][2] = rhs[2][2]/(1+4.0*C);
}
