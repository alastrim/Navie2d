#pragma once
#include "misc.h"

namespace fillers
{
double u1 (double t, double x, double y);
double u2 (double t, double x, double y);
double r (double t, double x, double y);
double p (double t, double x, double y);
double drdt (double t, double x, double y);
double dru1dx1 (double t, double x, double y);
double dru2dx2 (double , double x, double y);
double dru1dt (double t, double x, double y);
double dru2dt (double , double , double );
double dru1u1dx1 (double t, double x, double y);
double dru2u2dx2 (double t, double x, double y);
double dru2u1dx2 (double t, double x, double y);
double dru2u1dx1 (double t, double x, double y);
double dpdx1 (double t, double x, double y);
double dpdx2 (double t, double x, double y);
double f_first (double t, double x, double y);
double f_second (double t, double x, double y);
double f_third (double t, double x, double y);
}
