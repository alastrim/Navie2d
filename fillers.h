#ifndef FILLERS_H
#define FILLERS_H

#include <vector>
class discrete_function;

void fill_first (std::vector<double> &A, std::vector<double> &B, discrete_function &H, discrete_function &V, int n, double h, double tau, discrete_function &f, discrete_function &f_0);
void fill_second (std::vector<double> &A, std::vector<double> &B, discrete_function &H, discrete_function &V, int n, double h, double tau, discrete_function &f, discrete_function &f_0);

namespace fillers
{
double r (double t, double x);
double u (double t, double x);
double p (double t, double x);
double drdx (double t, double x);
double drdt (double t, double x);
double drudx (double t, double x);
double dudt (double t, double x);
double dudx (double t, double x);
double dpdx (double t, double x);
double ddudxx (double t, double x);
double f0 (double t, double x);
double f (double t, double x);
}

#endif // FILLERS_H
