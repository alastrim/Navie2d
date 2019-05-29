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
  assert (x > 0, "Sanity");

  return pow (x, y);
}
static double Pi = M_PI;
static double E = M_E;


namespace fillers
{
double u1 (double t, double x, double y)
{
  return sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (t);
}
double u2 (double t, double x, double y)
{
  return sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (-t);
}
double r (double t, double x, double y)
{
  return (cos (2.0 * M_PI * x) + 3.0/2.0) * (sin (2.0 * M_PI * y) + 3.0/2.0) * exp (t);
}
double p (double r)
{
  assert (r > 0, "Rho should be positive");
  return pow (r, GAMMA);
}

double dpdx (double t, double x, double y)
{
  return -2*Power(E,t)*GAMMA*Pi*Sin(2*Pi*x)*(1.5 + Sin(2*Pi*y))*Power(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)),-1 + GAMMA);
}


double f_1 (double t, double x, double y)
{
  auto du1dx = [](double t, double x, double y)
  { return 2. * Pi * cos (2. * Pi * x) * sin (2. * Pi * y) * exp (t); };
  auto du2dy = [](double t, double x, double y)
  { return 2. * Pi * sin (2. * Pi * x) * cos (2. * Pi * y) * exp (-t); };
  auto drdx = [](double t, double x, double y)
  { return - 2. * Pi * sin (2. * Pi * x) * (sin (2. * Pi * y) + 3./2.) * exp (t); };
  auto drdy = [](double t, double x, double y)
  { return 2. * Pi * (cos (2. * Pi * x) + 3./2.) * cos (2. * Pi * y) * exp (t); };
  auto drdt = [](double t, double x, double y)
  { return r (t, x, y); };

  return (drdt (t, x, y)) * 1.0
      + (r (t, x, y) * du1dx (t, x, y) + u1 (t, x, y) * drdx (t, x, y)) * 1.0
      + (r (t, x, y) * du2dy (t, x, y) + u2 (t, x, y) * drdy (t, x, y)) * 1.0;
}
double f_2 (double t, double x, double y)
{
  double res = (-2.*Power(E,t)*GAMMA*Pi*Sin(2.*Pi*x)*(1.5 + Sin(2.*Pi*y))*Power(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2.*Pi*y)),-1. + GAMMA) -
                MIU*((4.*Power(Pi,2.)*Cos(2.*Pi*x)*Cos(2.*Pi*y))/(3.*Power(E,t)) - (28.*Power(E,t)*Power(Pi,2.)*Sin(2*Pi*x)*Sin(2.*Pi*y))/3.) +
                Power(E,t)*(1.5 + Cos(2.*Pi*x))*(1.5 + Sin(2.*Pi*y))*(Power(E,t)*Sin(2.*Pi*x)*Sin(2.*Pi*y) + 2*Pi*Cos(2.*Pi*y)*Power(Sin(2.*Pi*x),2.)*Sin(2.*Pi*y) +
                   2.*Power(E,2.*t)*Pi*Cos(2.*Pi*x)*Sin(2.*Pi*x)*Power(Sin(2.*Pi*y),2.)))/(Power(E,t)*(1.5 + Cos(2.*Pi*x))*(1.5 + Sin(2.*Pi*y)));
  return res;
}
double f_3 (double t, double x, double y)
{
  double res = (2*Power(E,t)*GAMMA*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Power(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)),-1 + GAMMA) -
                MIU*((4*Power(E,t)*Power(Pi,2)*Cos(2*Pi*x)*Cos(2*Pi*y))/3. - (28*Power(Pi,2)*Sin(2*Pi*x)*Sin(2*Pi*y))/(3.*Power(E,t))) +
                Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y))*(-((Sin(2*Pi*x)*Sin(2*Pi*y))/Power(E,t)) + (2*Pi*Cos(2*Pi*y)*Power(Sin(2*Pi*x),2)*Sin(2*Pi*y))/Power(E,2*t) +
                   2*Pi*Cos(2*Pi*x)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)))/(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)));
  return res;
}
}
