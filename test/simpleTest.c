#include <cstdio>
#include <cmath>
#include <gsl/gsl_integration.h>
#include "rhoFun.h"

int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
  //    double sigma = *((double *) fdata); // we can pass σ via fdata argument
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    
    fval[0] = rho_p(x[0], x[1], R, b);
    return 0;
    /*

    double sum = 0;
    unsigned i;
    for (i = 0; i < ndim; ++i) sum += x[i] * x[i];
    // compute the output value: note that fdim should == 1 from below
    fval[0] = exp(-sigma * sum);
    return 0; // success
    */
}

int main() {
  /*
    double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5, val, err;
    hcubature(1, f, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
    printf("Computed integral = %0.10g +/- %g\n", val, err);
    return 0;
  */
  double R = 7.0; // radius of the nuclei R
  double b = 4.0; // impact parameter b
  double fdata[2] = {R,b};
  double xmin[2] = {-R-b/2,-R-b/2}, xmax[2] = {R+b/2,R+b/2};
  double val, err;
  hcubature(1, f, &fdata, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
  printf("Computed integral = %0.10g +/- %g\n", val, err);
  return 0;
}
