#include <math.h>
#include "fillers.h"
#include "misc.h"
#include "matvec.h"
#include "discrete_function.h"

namespace fillers
{
double r (double t, double x)
{
  if (GAP2)
    return 1;
  if (GAP)
    return ((x < 4.5 || x > 5.5) ? 1 : 2);
  return exp(t) * (cos(M_PI * x / 10.0) + 1.5);
}

double drdx (double t, double x)
{
  return - M_PI / 10 * exp(t) * sin (M_PI * x / 10.0);
}

double u (double t, double x)
{
  if (GAP2)
    return ((x < 4.5 || x > 5.5) ? 0 : 1);
  if (GAP)
    return 0;
  return cos (2.0 * M_PI * t) * sin (M_PI * (x / 10.0) * (x / 10.0));
}

double p (double t, double x)
{
  return pow (r (t, x), GAMMA);
}

double drdt (double t, double x)
{
  return exp(t) * (cos(M_PI * x / 10.0) + 1.5);
}

double drudx (double t, double x)
{
  return (exp (t) * M_PI * x * cos (2.0 * M_PI * t) * (1.5 + cos ((M_PI * x) / 10.0)) * cos ((M_PI * x * x) / 100.0)) / 50.0 - (exp (t) * M_PI * cos(2.0 * M_PI * t) * sin ((M_PI * x) / 10.0) * sin ((M_PI * x * x) / 100.0)) / 10.0;
}

double dudt (double t, double x)
{
  return -2.0 * M_PI * sin (2 * M_PI * t) * sin (M_PI * (x / 10.0) * (x / 10.0));
}

double dudx (double t, double x)
{
  return 1.0 / 50.0 * M_PI * x * cos (2 * M_PI * t) * cos (M_PI * (x / 10.0) * (x / 10.0));
}

double dpdx (double t, double x)
{
  return -(exp (t) * M_PI * GAMMA * pow ((exp (t) * (1.5 + cos ((M_PI * x) / 10.0))), (-1.0 + GAMMA)) * sin ((M_PI * x) / 10.0)) / 10.0;
}

double ddudxx (double t, double x)
{
  return (M_PI * cos (2.0 * M_PI * t) * cos ((M_PI * x * x) / 100.0))/50.0 - (M_PI * M_PI * x * x * cos (2.0 * M_PI * t) * sin ((M_PI * x * x) / 100.0)) / 2500.0;
}

double f0 (double t, double x)
{
  if (GAP)
    return 0;
  return drdt (t, x) + drudx (t, x);
}
double f (double t, double x)
{
  if (GAP)
    return 0;
  return (r (t, x) * dudt (t, x) + r (t, x) * u (t, x) * dudx (t, x) + dpdx (t, x) - MIU * ddudxx (t, x)) / r (t, x);
}
}

void fill_first (std::vector<double> &A, std::vector<double> &B,
                 discrete_function &H, discrete_function &V,
                 int n, double h, double tau, discrete_function &/*f*/,
                 discrete_function &f_0)
{
  std::vector<double> &H_cut = H.cut (n);
  std::vector<double> &V_cut = V.cut (n);
  unsigned int M = static_cast<unsigned int> (H_cut.size ());

  A = std::vector<double> (M * M, 0);
  B = std::vector<double> (M, 0);

  for (unsigned int m = 0; m < M; m++)
    {
      if (m > 0)
        A[m * M + m - 1] = -(V_cut[m] + fabs (V_cut[m])) / (2.0 * h);

      A[m * M + m] = 1.0 / tau + (V_cut[m + 1] + fabs (V_cut[m + 1])
          - V_cut[m] + fabs (V_cut[m])) / (2.0 * h);

      if (m < M - 1)
        A[m * M + m + 1] = (V_cut[m + 1] - fabs (V_cut[m + 1])) / (2.0 * h);

      B[m] = H_cut[m] * (1.0 / tau) + f_0.val (n, toi (m));
    }
}

void fill_second (std::vector<double> &A, std::vector<double> &B,
                  discrete_function &H, discrete_function &V,
                  int n, double h, double tau, discrete_function &f,
                  discrete_function &/*f_0*/)
{
  std::vector<double> &H_cut = H.cut (n + 1);
  std::vector<double> &V_cut = V.cut (n);
  unsigned int M = static_cast<unsigned int> (H_cut.size ());

  A = std::vector<double> ((M + 1) * (M + 1), 0);
  B = std::vector<double> (M + 1, 0);

  for (unsigned int m = 0; m < M + 1; m++)
    {
      double H_s_;
      if (m == 0 || m == M || !fuzzycmp ((H_s_ = (H_cut[m] + H_cut[m - 1]) / 2)))
        {
          A[m * (M + 1) + m] = 1;
          B[m] = 0;
          continue;
        }

      A[m * (M + 1) + m - 1] = H_s_ * -(V_cut[m] + fabs (V_cut[m]))
                               / (2 * h) - MIU / (h * h);

      A[m * (M + 1) + m] = H_s_ * (1.0 / tau + (V_cut[m] + fabs (V_cut[m]))
                                   / (2 * h) - (V_cut[m] - fabs (V_cut[m]))
                                   / (2 * h)) + (2 * MIU) / (h * h);

      A[m * (M + 1) + m + 1] = H_s_ * (V_cut[m] - fabs (V_cut[m]))
                               / (2 * h) - MIU / (h * h);

      B[m] = H_s_ * (V_cut[m] / tau - (GAMMA / (GAMMA - 1)) *
                     ((pow (H_cut[m], GAMMA - 1) - pow (H_cut[m - 1], GAMMA - 1))
                     / h) + f.val (n + 1, toi (m)));
    }
}
