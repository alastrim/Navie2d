﻿#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "misc.h"
#include "matvec.h"
#include "discrete_function.h"
#include "grid.h"


static double xmabs (double x)
{
  return x - fabs (x);
}

static double xpabs (double x)
{
  return x + fabs (x);
}

namespace fillers
{
void fill_first (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
  discrete_function &H = essential.m_tdfH.get_cut (k);
  discrete_function &V1 = essential.m_tdfV1.get_cut (k);
  discrete_function &V2 = essential.m_tdfV2.get_cut (k);
  grid_parameters g_parameters = essential.m_tdfV1.get_grid ()->get_parameters ();
  scale_parameters s_parameters = essential.m_tdfV1.get_scale ()->get_parameters ();

  unsigned int M1 = tou (g_parameters.m_x_step_count);
  unsigned int M2 = tou (g_parameters.m_y_step_count);
  unsigned int S = static_cast<unsigned int> (H.get_raw_vector ().size ());
  assert (S == M1 * M2, "Bad matrix size");

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  std::function<unsigned int (unsigned int i, unsigned int j)> vect_to_mat = [&] (unsigned int i, unsigned int j) { return i * S + j; };

  A = std::vector<double> (S * S, 0);
  B = std::vector<double> (S, 0);

  for (unsigned int m1 = 0; m1 < M1; m1++)
    {
      for (unsigned int m2 = 0; m2 < M2; m2++)
        {
          A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0/tau
                                               +(xpabs(V1.get_value({m1+1,m2})+V1.get_value({m1+1,m2+1}))
                                                  -xmabs(V1.get_value({m1,m2})+V1.get_value({m1,m2+1})))/2.0/h1
                                                 +(xpabs(V2.get_value({m1,m2+1})+V2.get_value({m1+1,m2+1}))
                                                   -xmabs(V2.get_value({m1,m2})+V2.get_value({m1+1,m2})))/2.0/h2;
          if (m2 < M2 - 1)
            A[vect_to_mat (m1*M2+m2,m1*M2+m2+1)] = xmabs((V2.get_value({m1,m2+1})+V2.get_value({m1+1,m2+1}))/2.0)/2.0/h2;
          if (m1 < M1 - 1)
            A[vect_to_mat (m1*M2+m2,(m1+1)*M2+m2)] = xmabs((V1.get_value({m1+1,m2})+V1.get_value({m1+1,m2+1}))/2.0)/2.0/h1;
          if (m1 > 0)
            A[vect_to_mat (m1*M2+m2,(m1-1)*M2+m2)] = -xpabs((V1.get_value({m1,m2})+V1.get_value({m1,m2+1}))/2.0)/2.0/h1;
          if (m2 > 0)
            A[vect_to_mat (m1*M2+m2,m1*M2+m2-1)] = -xpabs((V2.get_value({m1,m2})+V2.get_value({m1+1,m2}))/2.0)/2.0/h2;

          B[vect_to_mat (0,m1*M2+m2)] = 1.0/tau*H.get_value({m1,m2});
        }
    }
}

void fill_second (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
  discrete_function &H = essential.m_tdfH.get_cut (k + 1);
  discrete_function &V1 = essential.m_tdfV1.get_cut (k);
  discrete_function &V2 = essential.m_tdfV2.get_cut (k);
  grid_parameters g_parameters = essential.m_tdfV1.get_grid ()->get_parameters ();
  scale_parameters s_parameters = essential.m_tdfV1.get_scale ()->get_parameters ();

  unsigned int M1 = tou (g_parameters.m_x_step_count);
  unsigned int M2 = tou (g_parameters.m_y_step_count);
  unsigned int S = static_cast<unsigned int> (V1.get_raw_vector ().size ());
  assert (S == (M1 + 1) * (M2 + 1), "Bad matrix size");

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  std::function<unsigned int (unsigned int i, unsigned int j)> vect_to_mat = [&] (unsigned int i, unsigned int j) { return i * S + j; };

  A = std::vector<double> (S * S, 0);
  B = std::vector<double> (S, 0);

  for (unsigned int m1 = 0; m1 <= M1; m1++)
    {
      for (unsigned int m2 = 0; m2 <= M2; m2++)
        {
          if (m1 == 0 || m2 == 0 || m1 == M1 || m2 == M2)
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[vect_to_mat (0,m1*M2+m2)] = 0;
              continue;
            }

          double check = (H.get_value ({m1,m2})+H.get_value ({m1,m2-1})+H.get_value ({m1-1,m2})+H.get_value ({m1-1,m2-1}))/4;
          if (!fuzzycmp (check))
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[vect_to_mat (0,m1*M2+m2)] = 0;
              continue;
            }

          A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = check/tau
                                               +check*(+fabs(V1.get_value ({m1,m2}))/h1
                                                       +fabs(V2.get_value ({m1,m2}))/h2)
                                               +MIU*4.0/3.0*2/h1/h1 + MIU*1.0*2/h2/h2;
          if (m2 < M2 - 1)
            A[vect_to_mat (m1*M2+m2,m1*M2+m2+1)] = +check*xmabs(V2.get_value ({m1,m2}))/2.0/h2
                                                   -MIU*1.0/h2/h2;
          if (m1 < M1 - 1)
            A[vect_to_mat (m1*M2+m2,(m1+1)*M2+m2)] = +check*xmabs(V1.get_value ({m1,m2}))/2.0/h1
                                                     -MIU*4.0/3.0/h1/h1;
          if (m1 > 0)
            A[vect_to_mat (m1*M2+m2,(m1-1)*M2+m2)] = -check*xpabs(V1.get_value ({m1,m2}))/2.0/h1
                                                     -MIU*4.0/3.0/h1/h1;
          if (m2 > 0)
            A[vect_to_mat (m1*M2+m2,m1*M2+m2-1)] = -check*xpabs(V2.get_value ({m1,m2}))/2.0/h2
                                                   -MIU*1.0/h2/h2;

          B[vect_to_mat (0,m1*M2+m2)] = check/tau*V1.get_value ({m1,m2})
                                        -1.0/h1*(
                                                  +pow((H.get_value ({m1,m2})+H.get_value ({m1,m2-1}))/2,GAMMA)
                                                  -pow((H.get_value ({m1-1,m2})+H.get_value ({m1-1,m2-1}))/2,GAMMA))
                                        +MIU/3.0/4.0/h1/h2*(
                                                  +V2.get_value ({m1-1,m2-1})
                                                  -V2.get_value ({m1-1,m2+1})
                                                  -V2.get_value ({m1+1,m2-1})
                                                  +V2.get_value ({m1+1,m2+1}));
        }
    }
}

void fill_third (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
  discrete_function &H = essential.m_tdfH.get_cut (k + 1);
  discrete_function &V1 = essential.m_tdfV1.get_cut (k);
  discrete_function &V2 = essential.m_tdfV2.get_cut (k);
  grid_parameters g_parameters = essential.m_tdfV1.get_grid ()->get_parameters ();
  scale_parameters s_parameters = essential.m_tdfV1.get_scale ()->get_parameters ();

  unsigned int M1 = tou (g_parameters.m_x_step_count);
  unsigned int M2 = tou (g_parameters.m_y_step_count);
  unsigned int S = static_cast<unsigned int> (V1.get_raw_vector ().size ());
  assert (S == (M1 + 1) * (M2 + 1), "Bad matrix size");

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  std::function<unsigned int (unsigned int i, unsigned int j)> vect_to_mat = [&] (unsigned int i, unsigned int j) { return i * S + j; };

  A = std::vector<double> (S * S, 0);
  B = std::vector<double> (S, 0);

  for (unsigned int m1 = 0; m1 <= M1; m1++)
    {
      for (unsigned int m2 = 0; m2 <= M2; m2++)
        {
          if (m1 == 0 || m2 == 0 || m1 == M1 || m2 == M2)
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[vect_to_mat (0,m1*M2+m2)] = 0;
              continue;
            }

          double check = (H.get_value ({m1,m2})+H.get_value ({m1,m2-1})+H.get_value ({m1-1,m2})+H.get_value ({m1-1,m2-1}))/4;
          if (!fuzzycmp (check))
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[vect_to_mat (0,m1*M2+m2)] = 0;
              continue;
            }

          A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = check/tau
                                               +check*(+fabs(V1.get_value ({m1,m2}))/h1
                                                       +fabs(V2.get_value ({m1,m2}))/h2)
                                               +MIU*1.0*2/h1/h1 + MIU*4.0/3.0*2/h2/h2;
          if (m2 < M2 - 1)
            A[vect_to_mat (m1*M2+m2,m1*M2+m2+1)] = +check*xmabs(V2.get_value ({m1,m2}))/2.0/h2
                                                   -MIU*4.0/3.0/h2/h2;
          if (m1 < M1 - 1)
            A[vect_to_mat (m1*M2+m2,(m1+1)*M2+m2)] = +check*xmabs(V1.get_value ({m1,m2}))/2.0/h1
                                                     -MIU*1.0/h2/h2;
          if (m1 > 0)
            A[vect_to_mat (m1*M2+m2,(m1-1)*M2+m2)] = -check*xpabs(V1.get_value ({m1,m2}))/2.0/h1
                                                     -MIU*1.0/h1/h1;
          if (m2 > 0)
            A[vect_to_mat (m1*M2+m2,m1*M2+m2-1)] = -check*xpabs(V2.get_value ({m1,m2}))/2.0/h2
                                                   -MIU*4.0/3.0/h2/h2;

          B[vect_to_mat (0,m1*M2+m2)] = check/tau*V2.get_value ({m1,m2})
                                        -1.0/h2*(
                                                  +pow((H.get_value ({m1,m2})+H.get_value ({m1-1,m2}))/2,GAMMA)
                                                  -pow((H.get_value ({m1,m2-1})+H.get_value ({m1-1,m2-1}))/2,GAMMA))
                                        +MIU/3.0/4.0/h1/h2*(
                                                  +V1.get_value ({m1-1,m2-1})
                                                  -V1.get_value ({m1-1,m2+1})
                                                  -V1.get_value ({m1+1,m2-1})
                                                  +V1.get_value ({m1+1,m2+1}));
        }
    }
}
}
