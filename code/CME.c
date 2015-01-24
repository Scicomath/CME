#include <stdlib.h>
#include <stdio.h>
#include "head/eBpFun.h"
#include "head/eBsFun.h"
#include "head/eB.h"

const double alpha_EM = 1.0/137.0;
const double hbarc = 197.32696;
const size_t calls = 500000;
/*
double origin_eB(double t, double b, double Y0, double Z, double R, double a) {
  double x = 0.0, y = 0.0, z = 0.0;
  return 2 * eBs_p_y(x, y, z, t, b, Y0, Z, R) + 
    2 * eBp_p_y(x, y, z, t, b, Y0, Z, R, a);
}

double eB_x(double x, double y, double z, double t, double b, double Y0, double Z, double R, double a) {
  return eBs_p_x(x, y, z, t, b, Y0, Z, R) +
    eBs_m_x(x, y, z, t, b, Y0, Z, R) +
    eBp_p_x(x, y, z, t, b, Y0, Z, R, a) +
    eBp_m_x(x, y, z, t, b, Y0, Z, R, a);
}
*/
int main() {
  double t, b, Y0, Z, R, a;
  t = 0.5;
  b = 8;
  Y0 = 4.19;
  Z = 79.0;
  R = 7.0;
  a = 0.5;
  /*
  double eval1, eval2, eval3;
  // printf("-----------%f\n",eB_x(0, 1.0, 0, t,b,Y0,Z,R,a));
  
  for (t = 0.01;t<=3;t+=0.01) {
    b = 4.0;
    eval1 = origin_eB(t,b,Y0,Z,R,a);
    b = 8.0;
    eval2 = origin_eB(t,b,Y0,Z,R,a);
    b = 12.0;
    eval3 = origin_eB(t,b,Y0,Z,R,a);
    printf("%lg\t%lg\t%lg\n",eval1, eval2, eval3);
  }
  */
   
  double eval;            
  for (t = 0.01; t <= 1.0; t += 0.01) {
    printf("%lg\t",t);
  }
  printf("\n");
  for (Y0 = 0.01;Y0<=10;Y0+=0.1) {
    for (t = 0.01;t<=1.0;t+=0.01) {
      eval = origin_eB(t,b,Y0,Z,R,a);
      printf("%lg\t",eval);
    }
    printf("\n");
  }
  

  return 0;
}
