#include <stdlib.h>

#include "head/eBpFun.h"
#include "head/eBsFun.h"
#include "head/eB.h"

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
