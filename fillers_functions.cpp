#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "misc.h"
#include "matvec.h"
#include "discrete_function.h"
#include "grid.h"

static double Cos (double x) { return cos (x); }
static double Sin (double x) { return sin (x); }
static double Power (double x, double y)
{
  if (!fuzzycmp (y, 2))
    return x * x;
  if (!fuzzycmp (y, 3))
    return x * x * x;
  if (!fuzzycmp (x, 0))
    return 0;
  assert (x > 0, "X in power function should be positive");

  return pow (x, y);
}
static double E = M_E;


namespace fillers
{
double u1 (double t, double x, double y)
{
  return sin (x) * sin (y) * Power (E, t);
}
double u2 (double t, double x, double y)
{
  return sin (x) * sin (y) / Power (E, t);
}
double r (double t, double x, double y)
{
  return (cos (x) + 3.0/2.0) * (sin (y) + 3.0/2.0) * Power (E, t);
}
double p (double r)
{
  assert (r > 0, "Rho should be positive");
  return pow (r, GAMMA);
}

double f_1 (double t, double x, double y)
{
  return (1.5 + Cos(x))*Cos(y)*Sin(x)*Sin(y) + Power(E,t)*(1.5 + Cos(x))*(1.5 + Sin(y)) + (1.5 + Cos(x))*Cos(y)*Sin(x)*(1.5 + Sin(y)) + Power(E,2*t)*Cos(x)*(1.5 + Cos(x))*Sin(y)*(1.5 + Sin(y)) -
     Power(E,2*t)*Power(Sin(x),2)*Sin(y)*(1.5 + Sin(y));
}
double f_2 (double t, double x, double y)
{
  return (-(Power(E,t)*GAMMA*Sin(x)*(1.5 + Sin(y))*Power(Power(E,t)*(1.5 + Cos(x))*(1.5 + Sin(y)),-1 + GAMMA)) - MIU*((Cos(x)*Cos(y))/(3.*Power(E,t)) - (7*Power(E,t)*Sin(x)*Sin(y))/3.) +
          Power(E,t)*(1.5 + Cos(x))*(1.5 + Sin(y))*(Power(E,t)*Sin(x)*Sin(y) + Cos(y)*Power(Sin(x),2)*Sin(y) + Power(E,2*t)*Cos(x)*Sin(x)*Power(Sin(y),2)))/(Power(E,t)*(1.5 + Cos(x))*(1.5 + Sin(y)));
}
double f_3 (double t, double x, double y)
{
  return (Power(E,t)*GAMMA*(1.5 + Cos(x))*Cos(y)*Power(Power(E,t)*(1.5 + Cos(x))*(1.5 + Sin(y)),-1 + GAMMA) - MIU*((Power(E,t)*Cos(x)*Cos(y))/3. - (7*Sin(x)*Sin(y))/(3.*Power(E,t))) +
          Power(E,t)*(1.5 + Cos(x))*(1.5 + Sin(y))*(-((Sin(x)*Sin(y))/Power(E,t)) + (Cos(y)*Power(Sin(x),2)*Sin(y))/Power(E,2*t) + Cos(x)*Sin(x)*Power(Sin(y),2)))/(Power(E,t)*(1.5 + Cos(x))*(1.5 + Sin(y)));
}
}
