#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "misc.h"
#include "matvec.h"
#include "discrete_function.h"
#include "grid.h"

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
double drdt (double t, double x, double y)
{
  return exp (t) * (1.5 + cos (2.0 * M_PI * x)) * (1.5 + cos (2.0 * M_PI * y));
}
double dru1dx1 (double t, double x, double y)
{
  return 2.0 * exp (2.0 * t) * M_PI * cos (2.0 * M_PI * x) * (1.5 + cos (2.0 * M_PI * x)) * (1.5 + cos (2.0 * M_PI * y)) * sin (2.0 * M_PI * y) - 2.0 * exp (2.0 * t) * M_PI * (1.5 + cos (2.0 * M_PI * y)) * pow (sin (2.0 * M_PI * x), 2.0) * sin (2.0 * M_PI * y);
}
double dru2dx2 (double , double x, double y)
{
  return 2.0 * M_PI * (1.5 + cos (2.0 * M_PI * x)) * cos (2.0 * M_PI * y) * (1.5 + cos (2.0 * M_PI * y)) * sin (2.0 * M_PI * x) - 2.0 * M_PI * (1.5 + cos (2.0 * M_PI * x)) * sin (2.0 * M_PI * x) * pow (sin (2.0 * M_PI * y), 2.0);
}
double dru1dt (double t, double x, double y)
{
  return 2.0 * exp (2.0 * t) * (1.5 + cos (2.0 * M_PI * x)) * (1.5 + cos (2.0 * M_PI * y)) * sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y);
}
double dru2dt (double , double , double )
{
  return 0;
}
double dru1u1dx1 (double t, double x, double y)
{
  return 4.0 * exp (3.0 * t) * M_PI * cos (2.0 * M_PI * x) * (1.5 + cos (2.0 * M_PI * x)) * (1.5 + cos (2.0 * M_PI * y)) * sin (2.0 * M_PI * x) * pow (sin (2.0 * M_PI * y), 2.0) - 2.0 * exp (3.0 * t) * M_PI * (1.5 + cos (2.0 * M_PI * y)) * pow (sin (2.0 * M_PI * x), 3.0) * pow (sin (2.0 * M_PI * y), 2.0);
}
double dru2u2dx2 (double t, double x, double y)
{
  return 4.0 * exp (-t) * M_PI * (1.5 + cos (2.0 * M_PI * x)) * cos (2 * M_PI * y) * (1.5 + cos (2.0 * M_PI * y)) * pow (sin (2.0 * M_PI * x), 2.0) * sin (2.0 * M_PI * y) - 2.0 * exp (-t) * M_PI * (1.5 + cos (2 * M_PI * x)) * pow (sin (2.0 * M_PI * x), 2.0) * pow (sin (2.0 * M_PI * y), 3.0);
}
double dru2u1dx2 (double t, double x, double y)
{
  return 4.0 * exp (t) * M_PI * (1.5 + cos (2.0 * M_PI * x)) * cos (2.0 * M_PI * y) * (1.5 + cos (2.0 * M_PI * y)) * pow (sin (2.0 * M_PI * x), 2.0) * sin (2.0 * M_PI * y) - 2.0 * exp (t) * M_PI * (1.5 + cos (2.0 * M_PI * x)) * pow (sin (2.0 * M_PI * x), 2.0) * pow (sin (2.0 * M_PI * y), 3.0);
}
double dru2u1dx1 (double t, double x, double y)
{
  return 4.0 * exp (t) * M_PI * cos (2.0 * M_PI * x) * (1.5 + cos (2.0 * M_PI * x)) * (1.5 + cos (2.0 * M_PI * y)) * sin (2.0 * M_PI * x) * pow (sin (2.0 * M_PI * y), 2.0) - 2.0 * exp (t) * M_PI * (1.5 + cos (2.0 * M_PI * y)) * pow (sin (2.0 * M_PI * x), 3.0) * pow (sin (2.0 * M_PI * y), 2.0);
}
double dpdx1 (double t, double x, double y)
{
  return -2.0 * exp (t) * M_PI * GAMMA * (1.5 + cos (2.0 * M_PI * y)) * pow ((exp (t) * (1.5 + cos (2.0 * M_PI * x)) * (1.5 + cos (2.0 * M_PI * y))), (-1 + GAMMA)) * sin (2.0 * M_PI * x);
}
double dpdx2 (double t, double x, double y)
{
  return -2.0 * exp (t) * M_PI * GAMMA * (1.5 + cos (2.0 * M_PI * x)) * pow ((exp (t) * (1.5 + cos (2.0 * M_PI * x)) * (1.5 + cos (2.0 * M_PI * y))), (-1 + GAMMA)) *  sin (2.0 * M_PI * y);
}
double f_first (double t, double x, double y)
{
  return drdt (t, x, y) + dru1dx1 (t, x, y) + dru2dx2 (t, x, y);
}
double f_second (double t, double x, double y)
{
  return dru1dt (t, x, y) + dru1u1dx1 (t, x, y) + dru2u1dx2 (t, x, y) + dpdx1 (t, x, y);
}
double f_third (double t, double x, double y)
{
  return dru2dt (t, x, y) + dru2u1dx1 (t, x, y) + dru2u2dx2 (t, x, y) + dpdx2 (t, x, y);
}
}
