#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "misc.h"
#include "matvec.h"
#include "discrete_function.h"
#include "grid.h"

static double Cos (double x) { return cos (x); }
static double Sin (double x) { return sin (x); }
static double Power (double x, double y) { return pow (x, y); }
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






double f_1 (double t, double x, double y)
{
  double res = 2*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Sin(2*Pi*x)*Sin(2*Pi*y) +
               Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)) +
               2*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Sin(2*Pi*x)*(1.5 + Sin(2*Pi*y)) +
               2*Power(E,2*t)*Pi*Cos(2*Pi*x)*(1.5 + Cos(2*Pi*x))*Sin(2*Pi*y)*(1.5 + Sin(2*Pi*y)) -
               2*Power(E,2*t)*Pi*Power(Sin(2*Pi*x),2)*Sin(2*Pi*y)*(1.5 + Sin(2*Pi*y));

  return res;
}
double f_2 (double t, double x, double y)
{
  double res = (-2*Power(E,t)*GAMMA*Pi*Sin(2*Pi*x)*(1.5 + Sin(2*Pi*y))*
                Power(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)),-1 + GAMMA) -
               MIU*((4*Power(Pi,2)*Cos(2*Pi*x)*Cos(2*Pi*y))/(3.*Power(E,t)) -
                  (28*Power(E,t)*Power(Pi,2)*Sin(2*Pi*x)*Sin(2*Pi*y))/3.) +
               Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y))*
                (Power(E,t)*Sin(2*Pi*x)*Sin(2*Pi*y) + 2*Pi*Cos(2*Pi*y)*Power(Sin(2*Pi*x),2)*Sin(2*Pi*y) +
                  2*Power(E,2*t)*Pi*Cos(2*Pi*x)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)))/
             (Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)));
  return res;
}
double f_3 (double t, double x, double y)
{
  double res = (2*Power(E,t)*GAMMA*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*
                Power(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)),-1 + GAMMA) -
               MIU*((4*Power(E,t)*Power(Pi,2)*Cos(2*Pi*x)*Cos(2*Pi*y))/3. -
                  (28*Power(Pi,2)*Sin(2*Pi*x)*Sin(2*Pi*y))/(3.*Power(E,t))) +
               Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y))*
                (-((Sin(2*Pi*x)*Sin(2*Pi*y))/Power(E,t)) +
                  (2*Pi*Cos(2*Pi*y)*Power(Sin(2*Pi*x),2)*Sin(2*Pi*y))/Power(E,2*t) +
                  2*Pi*Cos(2*Pi*x)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)))/(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)));
  return res;
}
}
