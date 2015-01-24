#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "head/rhoFun.h"
#include "head/eBsFun.h"
#include "head/constant.h"


const size_t warmCalls = 100000;

double eBs_p_x(double x, double y, double z, double t, double b, double Y0, double Z, double R) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[8] = {R, b, x, y, z, tau, eta, Y0};
  double xmin[2] = {-R-b/2,-R}, xmax[2] = {R+b/2,R};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBs_p_x_int, 2, &fdata};

  // size_t calls = 500000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(2);
  gsl_monte_plain_integrate(&G, xmin, xmax, 2, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */

  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 2, warmCalls, r, s,
                               &intval, &err);

    //    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 2, calls/5, r, s,
                                   &intval, &err);
	//       printf ("result = % .6f sigma = % .6f "
	//      "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);

  /*
  hcubature(1, eBs_p_x_int, &fdata, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
  // printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = hbarc*hbarc * Z * alpha_EM * sinh(Y0-eta) * intval;
  // printf("Computed constant = %0.10g \n",hbarc*hbarc * Z * alpha_EM * sinh(Y0-eta) );
  return val;
}

double eBs_p_x_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] <= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] >= R*R)) {
	eval = rho_p(X_p[0], X_p[1], R, b) * (X_p[1] - y) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(Y0-eta),2.0) ,1.5);
      } else {
      eval = 0.0;
    }
    return eval;
}


double eBs_p_y(double x, double y, double z, double t, double b, double Y0, double Z, double R) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[8] = {R, b, x, y, z, tau, eta, Y0};
  double xmin[2] = {-R-b/2,-R}, xmax[2] = {R+b/2,R};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBs_p_y_int, 2, &fdata};

  // size_t calls = 500000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(2);
  gsl_monte_plain_integrate(&G, xmin, xmax, 2, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */

  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 2, warmCalls, r, s,
                               &intval, &err);

    // printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 2, calls/5, r, s,
                                   &intval, &err);
        // printf ("result = % .6f sigma = % .6f "
	//         "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);

  /*
  hcubature(1, eBs_p_y_int, &fdata, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
  // printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = hbarc*hbarc * Z * alpha_EM * sinh(Y0-eta) * intval;
  // printf("Computed constant = %0.10g \n",hbarc*hbarc * Z * alpha_EM * sinh(Y0-eta) );
  return val;
}

double eBs_p_y_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] <= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] >= R*R)) {
	eval = rho_p(X_p[0], X_p[1], R, b) * (x - X_p[0]) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(Y0-eta),2.0) ,1.5);
      } else {
      eval = 0;
    }

    return eval;
}


/* calculate eBs_m */

double eBs_m_x(double x, double y, double z, double t, double b, double Y0, double Z, double R) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[8] = {R, b, x, y, z, tau, eta, Y0};
  double xmin[2] = {-R-b/2,-R}, xmax[2] = {R+b/2,R};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBs_m_x_int, 2, &fdata};

  // size_t calls = 500000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(2);
  gsl_monte_plain_integrate(&G, xmin, xmax, 2, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */

  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 2, warmCalls, r, s,
                               &intval, &err);

    // printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 2, calls/5, r, s,
                                   &intval, &err);
	//  printf ("result = % .6f sigma = % .6f "
	//       "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);

  /*
  hcubature(1, eBs_m_x_int, &fdata, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
  // printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = - hbarc*hbarc * Z * alpha_EM * sinh(Y0+eta) * intval;
  // printf("Computed constant = %0.10g \n",- hbarc*hbarc * Z * alpha_EM * sinh(Y0+eta) );
  return val;
}

double eBs_m_x_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] >= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] <= R*R)) {
	eval = rho_m(X_p[0], X_p[1], R, b) * (X_p[1] - y) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(Y0+eta),2.0) ,1.5);
      } else {
      eval = 0;
    }

    return eval;
}


double eBs_m_y(double x, double y, double z, double t, double b, double Y0, double Z, double R) {
  double eta = 0.5 * log((t+z)/(t-z));
  double tau = pow(t*t - z*z,0.5);
  double fdata[8] = {R, b, x, y, z, tau, eta, Y0};
  double xmin[2] = {-R-b/2,-R}, xmax[2] = {R+b/2,R};
  double intval, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&eBs_m_y_int, 2, &fdata};

  // size_t calls = 500000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  /*
  gsl_monte_plain_state *s = gsl_monte_plain_alloc(2);
  gsl_monte_plain_integrate(&G, xmin, xmax, 2, calls, r, s, &intval, &err);
  gsl_monte_plain_free(s);
  */
  // monte vegas
  
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&G, xmin, xmax, 2, warmCalls, r, s,
                               &intval, &err);

    // printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xmin, xmax, 2, calls/5, r, s,
                                   &intval, &err);
        //printf ("result = % .6f sigma = % .6f "
        //        "chisq/dof = %.1f\n", intval, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s);
    gsl_rng_free(r);
  /*
  hcubature(1, eBs_m_y_int, &fdata, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &intval, &err);
  */
  // printf("Computed integral = %0.10g +/- %g\n", intval, err);
  double val = - hbarc*hbarc * Z * alpha_EM * sinh(Y0+eta) * intval;
  // printf("Computed constant = %0.10g \n",- hbarc*hbarc * Z * alpha_EM * sinh(Y0+eta) );
  return val;
}

double eBs_m_y_int(const double *X_p, size_t dim, void *fdata) {
    double R = ((double *) fdata)[0];
    double b = ((double *) fdata)[1];
    double x = ((double *) fdata)[2];
    double y = ((double *) fdata)[3];
    double z = ((double *) fdata)[4];
    double tau = ((double *) fdata)[5];
    double eta = ((double *) fdata)[6];
    double Y0 = ((double *) fdata)[7];
    double eval;

    if ((pow(X_p[0]+b/2,2.0) + X_p[1]*X_p[1] >= R*R)&&(pow(X_p[0]-b/2,2.0) + X_p[1]*X_p[1] <= R*R)) {
	eval = rho_m(X_p[0], X_p[1], R, b) * (x - X_p[0]) / pow( pow(X_p[0]-x,2.0) + pow(X_p[1]-y,2.0) + tau*tau * pow(sinh(Y0+eta),2.0) ,1.5);
      } else {
      eval = 0;
    }

    return eval;
}

