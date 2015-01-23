#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "head/rhoFun.h"
#include "head/eBpFun.h"
#include "head/constant.h"

// const double alpha_EM = 1.0/137.0;
// const double hbarc = 197.32696;
// const size_t calls = 2000000;
/*
int main() {
  double x = 0.0, y = 0.0, z = 0.0, t = 1.0, b = 12, Y0=4.2, Z=79.0, R=7.0, a=0.5;
  double val = eBp_p_y(x, y, z, t, b, Y0, Z, R, a);
  
  printf("Computed result = %0.10g \n", val);
  return 0;
}
*/
double eBp_p_x(double x, double y, double z, double t, double b, double Y0, double Z, double R, double a) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[9] = {R, b, x, y, z, tau, eta, Y0, a};
  double xmin[3] = {-R+b/2,-R,-Y0}, xmax[3] = {R-b/2,R,Y0};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBp_p_x_int, 3, &fdata};

  // size_t calls = 1000000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  // monte plain
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
  gsl_monte_plain_integrate(&G, xmin, xmax, 3, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */
  // monte miser
  /*
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G, xmin, xmax, 3, calls, r, s,
                               &intval, &err);
    gsl_monte_miser_free (s);
  */
  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 3, 10000, r, s,
                               &intval, &err);

    //    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 3, calls/5, r, s,
                                   &intval, &err);
	//  printf ("result = % .6f sigma = % .6f "
	//      "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);
  /*
  hcubature(1, eBp_p_x_int, &fdata, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
  //printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = hbarc*hbarc * Z * alpha_EM * intval;
  //printf("Computed constant = %0.10g \n",hbarc*hbarc * Z * alpha_EM );
  return val;
}

double eBp_p_x_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double a = ((double *) fdata)[8];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] <= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] <= R*R)) {
      eval = sinh(X_p[2]-eta) * (a/(2*sinh(a*Y0))*exp(a*X_p[2])) * rho_p(X_p[0], X_p[1], R, b) * (X_p[1] - y) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(X_p[2]-eta),2.0) ,1.5);
      //   printf("fval = %f\n", fval[0]);
      } else {
      eval = 0;
    }

    return eval;
}


double eBp_p_y(double x, double y, double z, double t, double b, double Y0, double Z, double R, double a) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[9] = {R, b, x, y, z, tau, eta, Y0, a};
  double xmin[3] = {-R-b/2,-R, -Y0}, xmax[3] = {R+b/2,R,Y0};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBp_p_y_int, 3, &fdata};

  // size_t calls = 1000000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
  gsl_monte_plain_integrate(&G, xmin, xmax, 3, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */

  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 3, 10000, r, s,
                               &intval, &err);

    // printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 3, calls/5, r, s,
                                   &intval, &err);
        // printf ("result = % .6f sigma = % .6f "
	//               "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);

  /*
  hcubature(1, eBp_p_y_int, &fdata, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
  // printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = hbarc*hbarc * Z * alpha_EM * intval;
// printf("Computed constant = %0.10g \n",hbarc*hbarc * Z * alpha_EM );
  return val;
}

double eBp_p_y_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double a = ((double *) fdata)[8];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] <= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] <= R*R)) {
      eval = sinh(X_p[2]-eta) * (a/(2*sinh(a*Y0))*exp(a*X_p[2])) * rho_p(X_p[0], X_p[1], R, b) * (x - X_p[0]) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(X_p[2]-eta),2.0) ,1.5);
      } else {
      eval = 0;
    }

    return eval;
}

/* calculate eBp_m */

double eBp_m_x(double x, double y, double z, double t, double b, double Y0, double Z, double R, double a) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[9] = {R, b, x, y, z, tau, eta, Y0, a};
  double xmin[3] = {-R-b/2,-R, -Y0}, xmax[3] = {R+b/2,R, Y0};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBp_m_x_int, 3, &fdata};

  // size_t calls = 1000000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
  gsl_monte_plain_integrate(&G, xmin, xmax, 3, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */

  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 3, 10000, r, s,
                               &intval, &err);

    // printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 3, calls/5, r, s,
                                   &intval, &err);
	//  printf ("result = % .6f sigma = % .6f "
	//      "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);

  /*
  hcubature(1, eBp_m_x_int, &fdata, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
    // printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = - hbarc*hbarc * Z * alpha_EM * intval;
  // printf("Computed constant = %0.10g \n",- hbarc*hbarc * Z * alpha_EM  );
  return val;
}

double eBp_m_x_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double a = ((double *) fdata)[8];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] <= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] <= R*R)) {
      eval = sinh(X_p[2]+eta) * (a/(2*sinh(a*Y0))*exp(a*X_p[2])) * rho_m(X_p[0], X_p[1], R, b) * (X_p[1] - y) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(X_p[2]+eta),2.0) ,1.5);
      } else {
      eval = 0;
    }

    return eval;
}


double eBp_m_y(double x, double y, double z, double t, double b, double Y0, double Z, double R, double a) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[9] = {R, b, x, y, z, tau, eta, Y0, a};
  double xmin[3] = {-R-b/2,-R,-Y0}, xmax[3] = {R+b/2,R,Y0};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBp_m_y_int, 3, &fdata};

  // size_t calls = 1000000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(3);
  gsl_monte_plain_integrate(&G, xmin, xmax, 3, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */

  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 3, 10000, r, s,
                               &intval, &err);

    //printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 3, calls/5, r, s,
                                   &intval, &err);
	//  printf ("result = % .6f sigma = % .6f "
	//      "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);

  /*
  hcubature(1, eBp_m_y_int, &fdata, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
    //  printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = - hbarc*hbarc * Z * alpha_EM * intval;
  //printf("Computed constant = %0.10g \n",- hbarc*hbarc * Z * alpha_EM );
  return val;
}

double eBp_m_y_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double a = ((double *) fdata)[8];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] <= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] <= R*R)) {
	eval = sinh(X_p[2]+eta) * (a/(2*sinh(a*Y0))*exp(a*X_p[2])) * rho_m(X_p[0], X_p[1], R, b) * (x - X_p[0]) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(X_p[2]+eta),2.0) ,1.5);
      } else {
      eval = 0;
    }

    return eval;
}
