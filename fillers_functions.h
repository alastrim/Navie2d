#pragma once
#include "misc.h"

namespace fillers
{
double u1 (double t, double x, double y);
double u2 (double t, double x, double y);
double r (double t, double x, double y);
double p (double t, double x, double y);
double ddu1dx1dx1 (double t, double x, double y);
double ddu1dx1dx2 (double t, double x, double y);
double ddu1dx2dx2 (double t, double x, double y);
double ddu2dx1dx1 (double t, double x, double y);
double ddu2dx1dx2 (double t, double x, double y);
double ddu2dx2dx2 (double t, double x, double y);
double drdt (double t, double x, double y);
double du1dx1 (double t, double x, double y);
double du2dx2 (double , double x, double y);
double du1dx2 (double t, double x, double y);
double du2dx1 (double , double x, double y);
double dru1dx1 (double t, double x, double y);
double dru2dx2 (double , double x, double y);
double du1dt (double t, double x, double y);
double du2dt (double , double , double );
double dpdx1 (double t, double x, double y);
double dpdx2 (double t, double x, double y);
double f_first (double t, double x, double y);
double f_second (double t, double x, double y);
double f_third (double t, double x, double y);
}
