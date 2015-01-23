#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include "head/eB.h"

const double alpha_EM = 1.0/137.0;
const double hbarc = 197.32696;
const size_t calls = 500000;

double fn1 (double Y0, void * params)
{
  double t = ((double *) params)[0];
  double b = ((double *) params)[1];
  double Z = ((double *) params)[2];
  double R = ((double *) params)[3];
  double a = ((double *) params)[4];

  return -origin_eB(t,b,Y0,Z,R,a);
}

double fminimizer(double a, double b, double m, double params[])
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;

  gsl_function F;
  
  
  F.function = &fn1;
  F.params = params;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc(T);
  gsl_min_fminimizer_set(s, &F, m, a, b);

  //printf("using %s method\n", gsl_min_fminimizer_name(s));

  //printf ("%5s [%9s, %9s] %9s %9s\n",
  //	    "iter", "lower", "upper", "min",
  //    "err(est)");

  //printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
  //	    iter, a, b,
  //	    m, b - a);

    do {
      iter++;
      status = gsl_min_fminimizer_iterate(s);

      m = gsl_min_fminimizer_x_minimum(s);
      a = gsl_min_fminimizer_x_lower(s);
      b = gsl_min_fminimizer_x_upper(s);

      status = gsl_min_test_interval(a,b,0.001,0.0);

      //if (status == GSL_SUCCESS)
      //printf("Converged:\n");

      //printf ("%5d [%.7f, %.7f] "
      //        "%.7f %.7f\n",
      //      iter, a, b,
      //      m, b - a);
    } while(status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free(s);

  return m;
}

int main(void)
{
  double m = 5.0;
  double a = 0.10, b = 10.0;
  double params[5] = {0.01,8.0,79.0,7.0,0.5};//t,b,Z,R,a
  double fmin;

  for(params[0]=0.01;params[0]<=1.0;params[0]+=0.01) {
    fmin = fminimizer(a, b, m, params);
    printf("%lg\t%lg\n",fmin,
	   origin_eB(params[0],params[1],fmin,params[2],params[3],params[4]));
  }
  
  return 0;
}




