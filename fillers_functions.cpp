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
double p (double t, double x, double y)
{
  return pow (r (t, x, y), GAMMA);
}
double ddu1dx1dx1 (double t, double x, double y)
{
  return - 4 * M_PI * M_PI * sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (t);
}
double ddu1dx1dx2 (double t, double x, double y)
{
  return 4 * M_PI * M_PI * cos (2.0 * M_PI * x) * cos (2.0 * M_PI * y) * exp (t);
}
double ddu1dx2dx2 (double t, double x, double y)
{
  return - 4 * M_PI * M_PI * sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (t);
}
double ddu2dx1dx1 (double t, double x, double y)
{
  return - 4 * M_PI * M_PI * sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (-t);
}
double ddu2dx1dx2 (double t, double x, double y)
{
  return 4 * M_PI * M_PI * cos (2.0 * M_PI * x) * cos (2.0 * M_PI * y) * exp (-t);
}
double ddu2dx2dx2 (double t, double x, double y)
{
  return - 4 * M_PI * M_PI * sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (-t);
}
double du1dx1 (double t, double x, double y)
{
  return 2 * M_PI * cos (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (t);
}
double du2dx2 (double t, double x, double y)
{
  return 2 * M_PI * sin (2.0 * M_PI * x) * cos (2.0 * M_PI * y) * exp (-t);
}
double du1dx2 (double t, double x, double y)
{
  return 2 * M_PI * sin (2.0 * M_PI * x) * cos (2.0 * M_PI * y) * exp (t);
}
double du2dx1 (double t, double x, double y)
{
  return 2 * M_PI * cos (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * exp (-t);
}
double drdt (double t, double x, double y)
{
  return r (t, x, y);
}
double dru1dx1 (double t, double x, double y)
{
  return 2.0 * exp (2.0 * t) * M_PI * (1.5 + sin (2.0 * M_PI * y)) * sin (2.0 * M_PI * y) *
      (- sin (2.0 * M_PI * x) * sin (2.0 * M_PI * x) + (cos (2.0 * M_PI * x) + 1.5) * sin (2.0 * M_PI * x));
}
double dru2dx2 (double , double x, double y)
{
  return 2.0 * M_PI * (1.5 + cos (2.0 * M_PI * x)) * sin (2.0 * M_PI * x) *
      (+ sin (2.0 * M_PI * y) * cos (2.0 * M_PI * y) + (sin (2.0 * M_PI * y) + 1.5) * sin (2.0 * M_PI * y));
}
double du1dt (double t, double x, double y)
{
  return u1 (t, x, y);
}
double du2dt (double t, double x, double y)
{
  return -u2 (t, x, y);
}
double dpdx1 (double t, double x, double y)
{
  return -2.0 * exp (t) * M_PI * GAMMA * (1.5 + sin (2.0 * M_PI * y)) * pow ((exp (t) * (1.5 + cos (2.0 * M_PI * x))
         * (1.5 + sin (2.0 * M_PI * y))), (-1 + GAMMA)) * sin (2.0 * M_PI * x);
}
double dpdx2 (double t, double x, double y)
{
  return 2.0 * exp (t) * M_PI * GAMMA * (1.5 + cos (2.0 * M_PI * x)) * pow ((exp (t) * (1.5 + cos (2.0 * M_PI * x))
         * (1.5 + sin (2.0 * M_PI * y))), (-1 + GAMMA)) *  cos (2.0 * M_PI * y);
}
double f_first (double t, double x, double y)
{
  double a = drdt (t, x, y);
  double b = dru1dx1 (t, x, y);
  double c = dru2dx2 (t, x, y);
  double res = a + b + c;
  return res;
}
double f_second (double t, double x, double y)
{
  double res = (r (t, x, y) * (du1dt (t, x, y) + u1 (t, x, y) * du1dx1 (t, x, y) + u2 (t, x, y) * du1dx2 (t, x, y))
      + dpdx1 (t, x, y)
      - MIU * (4.0/3.0 * ddu1dx1dx1 (t, x, y) + ddu1dx2dx2 (t, x, y) + 1.0/3.0 * ddu2dx1dx2 (t, x, y))) / r (t, x, y);
  return res;
}
double f_third (double t, double x, double y)
{
  double res =  (r (t, x, y) * (du2dt (t, x, y) + u1 (t, x, y) * du2dx1 (t, x, y) + u2 (t, x, y) * du2dx2 (t, x, y))
      + dpdx2 (t, x, y)
      - MIU * (4.0/3.0 * ddu2dx2dx2 (t, x, y) + ddu2dx1dx1 (t, x, y) + 1.0/3.0 * ddu1dx1dx2 (t, x, y))) / r (t, x, y);
  return res;
}
double f_1 (double t, double x, double y)
{
  double res = 2*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Sin(2*Pi*x)*Sin(2*Pi*y) +
               Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)) +
               2*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Sin(2*Pi*x)*(1.5 + Sin(2*Pi*y)) +
               2*Power(E,2*t)*Pi*Cos(2*Pi*x)*(1.5 + Cos(2*Pi*x))*Sin(2*Pi*y)*
                (1.5 + Sin(2*Pi*y)) - 2*Power(E,2*t)*Pi*Power(Sin(2*Pi*x),2)*Sin(2*Pi*y)*
                (1.5 + Sin(2*Pi*y));
//  printf ("diff1 = %e\n", res - f_first (t,x,y));
  return res;
}
double f_2 (double t, double x, double y)
{
  double res = (2*Power(E,t)*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Power(Sin(2*Pi*x),2)*
                Power(Sin(2*Pi*y),2) + 2*Power(E,2*t)*(1.5 + Cos(2*Pi*x))*Sin(2*Pi*x)*
                Sin(2*Pi*y)*(1.5 + Sin(2*Pi*y)) +
               4*Power(E,t)*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Power(Sin(2*Pi*x),2)*
                Sin(2*Pi*y)*(1.5 + Sin(2*Pi*y)) +
               4*Power(E,3*t)*Pi*Cos(2*Pi*x)*(1.5 + Cos(2*Pi*x))*Sin(2*Pi*x)*
                Power(Sin(2*Pi*y),2)*(1.5 + Sin(2*Pi*y)) -
               2*Power(E,3*t)*Pi*Power(Sin(2*Pi*x),3)*Power(Sin(2*Pi*y),2)*
                (1.5 + Sin(2*Pi*y)) - 2*Power(E,t)*GAMMA*Pi*Sin(2*Pi*x)*(1.5 + Sin(2*Pi*y))*
                Power(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)),-1 + GAMMA) -
               MIU*((4*Power(Pi,2)*Cos(2*Pi*x)*Cos(2*Pi*y))/(3.*Power(E,t)) -
                  (28*Power(E,t)*Power(Pi,2)*Sin(2*Pi*x)*Sin(2*Pi*y))/3.))/
             (Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)));
//  printf ("diff2 = %e\n", res - f_second (t,x,y));
  return res;
}
double f_3 (double t, double x, double y)
{
  double res = ((2*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Power(Sin(2*Pi*x),2)*Power(Sin(2*Pi*y),2))/
                Power(E,t) + (4*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*Power(Sin(2*Pi*x),2)*
                  Sin(2*Pi*y)*(1.5 + Sin(2*Pi*y)))/Power(E,t) +
               4*Power(E,t)*Pi*Cos(2*Pi*x)*(1.5 + Cos(2*Pi*x))*Sin(2*Pi*x)*
                Power(Sin(2*Pi*y),2)*(1.5 + Sin(2*Pi*y)) -
               2*Power(E,t)*Pi*Power(Sin(2*Pi*x),3)*Power(Sin(2*Pi*y),2)*(1.5 + Sin(2*Pi*y)) +
               2*Power(E,t)*GAMMA*Pi*(1.5 + Cos(2*Pi*x))*Cos(2*Pi*y)*
                Power(Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)),-1 + GAMMA) -
               MIU*((4*Power(E,t)*Power(Pi,2)*Cos(2*Pi*x)*Cos(2*Pi*y))/3. -
                  (28*Power(Pi,2)*Sin(2*Pi*x)*Sin(2*Pi*y))/(3.*Power(E,t))))/
             (Power(E,t)*(1.5 + Cos(2*Pi*x))*(1.5 + Sin(2*Pi*y)));
//  printf ("diff3 = %e\n", res - f_third (t,x,y));
  return res;
}
}
