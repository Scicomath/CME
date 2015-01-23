#include <math.h>
#include "head/rhoFun.h"

/* Number densities of the two nuclei */
const double pi = 3.14159265;

double rho_p(double x_p, double y_p, double R, double b ) {
  double rho = 0;
  double temp = R*R - (pow(x_p + b/2,2.0) + y_p*y_p );
  if (temp >= 0)
    rho = 2/(4*pi*pow(R,3.0)/3) * sqrt(temp);

  return rho;
}
double rho_m(double x_p, double y_p, double R, double b ) {
  double rho = 0;
  double temp = R*R - (pow(x_p - b/2,2.0) + y_p*y_p );
  if (temp >= 0)
    rho = 2/(4*pi*pow(R,3.0)/3) * sqrt(temp);

  return rho;
}
