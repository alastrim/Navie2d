#include "fillers_misc.h"
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
  const grid *gr = essential.m_tdfH.get_grid ();
  const scale *sc = essential.m_tdfH.get_scale ();
  grid_parameters g_parameters = gr->get_parameters ();
  scale_parameters s_parameters = sc->get_parameters ();

  unsigned int M1 = tou (g_parameters.m_x_step_count) + 1;
  unsigned int M2 = tou (g_parameters.m_y_step_count) + 1;
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
          double t = sc->get_time (k + 1), x = gr->get_point ({m1,m2}).first, y = gr->get_point ({m1,m2}).second;

          A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.
                                               +tau*(xpabs(V1.tilda(m1+1,m2))-xmabs(V1.tilda(m1,m2)))/2.0/h1
                                               +tau*(xpabs(V2.tilda(m1,m2+1))-xmabs(V2.tilda(m1,m2)))/2.0/h2;

          if (m2<M2-1)
            A[vect_to_mat(m1*M2+m2,m1*M2+m2+1)]=tau*xmabs(V2.tilda (m1,m2+1))/2.0/h2;
          if (m1<M1-1)
            A[vect_to_mat(m1*M2+m2,(m1+1)*M2+m2)]=tau*xmabs(V1.tilda (m1+1,m2))/2.0/h1;
          if (m1>0)
            A[vect_to_mat(m1*M2+m2,(m1-1)*M2+m2)]=-tau*xpabs(V1.tilda (m1,m2))/2.0/h1;
          if (m2>0)
            A[vect_to_mat(m1*M2+m2,m1*M2+m2-1)]=-tau*xpabs(V2.tilda (m1,m2))/2.0/h2;

          B[m1*M2+m2]=H.val(m1,m2)+tau*f_1(t,x,y);
        }
    }
}

void fill_second (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
  discrete_function &H = essential.m_tdfH.get_cut (k + 1);
  discrete_function &V1 = essential.m_tdfV1.get_cut (k);
  discrete_function &V2 = essential.m_tdfV2.get_cut (k);
  const grid *gr = essential.m_tdfV1.get_grid ();
  const scale *sc = essential.m_tdfV1.get_scale ();
  grid_parameters g_parameters = gr->get_parameters ();
  scale_parameters s_parameters = sc->get_parameters ();

  unsigned int M1 = tou (g_parameters.m_x_step_count) + 1;
  unsigned int M2 = tou (g_parameters.m_y_step_count) + 1;
  unsigned int S = static_cast<unsigned int> (V1.get_raw_vector ().size ());
  assert (S == (M1) * (M2), "Bad matrix size");

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
          double t = sc->get_time (k + 1), x = gr->get_point ({m1,m2}).first, y = gr->get_point ({m1,m2}).second;

          if (m1 == 0 || m2 == 0 || m1 == (M1 - 1) || m2 == (M2 - 1))
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[m1*M2+m2] = 0.0;
              continue;
            }

          double check = (H.val(m1,m2)+H.val(m1,m2-1)+H.val(m1-1,m2)+H.val(m1-1,m2-1))/4.0;
          if (!fuzzycmp (check))
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[m1*M2+m2] = 0.0;
              continue;
            }

          A[vect_to_mat(m1*M2+m2,m1*M2+m2)]=check*(1.+tau/h1*fabs(V1.val(m1,m2))+tau/h2*fabs(V2.val(m1,m2)))
                                            +tau*MIU*(8./3./h1/h1+2./h2/h2);

          A[vect_to_mat(m1*M2+m2,(m1-1)*M2+m2)]=-tau/2./h1*xpabs(V1.val(m1,m2))*check
                                                -4./3.*tau*MIU/h1/h1;

          A[vect_to_mat(m1*M2+m2,(m1+1)*M2+m2)]=tau/2./h1*xmabs(V1.val(m1,m2))*check
                                                -4./3.*tau*MIU/h1/h1;

          A[vect_to_mat(m1*M2+m2,m1*M2+m2-1)]=-tau/2./h2*xpabs(V2.val(m1,m2))*check
                                              -tau*MIU/h2/h2;

          A[vect_to_mat(m1*M2+m2,m1*M2+m2+1)]=tau/2./h2*xmabs(V2.val(m1,m2))*check
                                              -tau*MIU/h2/h2;

          B[m1*M2+m2]=check*V1.val(m1,m2)-tau/h1*(p(H.left(m1,m2))-p(H.left(m1-1,m2)))
                      +tau*MIU/12./h1/h2*(V2.val(m1+1,m2+1)-V2.val(m1+1,m2-1)-V2.val(m1-1,m2+1)+V2.val(m1-1,m2-1))
                      +tau*f_2(t,x,y)*check;
        }
    }
}

void fill_third (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
  discrete_function &H = essential.m_tdfH.get_cut (k + 1);
  discrete_function &V1 = essential.m_tdfV1.get_cut (k);
  discrete_function &V2 = essential.m_tdfV2.get_cut (k);
  const grid *gr = essential.m_tdfV1.get_grid ();
  const scale *sc = essential.m_tdfV1.get_scale ();
  grid_parameters g_parameters = gr->get_parameters ();
  scale_parameters s_parameters = sc->get_parameters ();

  unsigned int M1 = tou (g_parameters.m_x_step_count) + 1;
  unsigned int M2 = tou (g_parameters.m_y_step_count) + 1;
  unsigned int S = static_cast<unsigned int> (V2.get_raw_vector ().size ());
  assert (S == (M1) * (M2), "Bad matrix size");

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
          double t = sc->get_time (k + 1), x = gr->get_point ({m1,m2}).first, y = gr->get_point ({m1,m2}).second;

          if (m1 == 0 || m2 == 0 || m1 == (M1 - 1) || m2 == (M2 - 1))
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[m1*M2+m2] = 0.0;
              continue;
            }

          double check = (H.val(m1,m2)+H.val(m1,m2-1)+H.val(m1-1,m2)+H.val(m1-1,m2-1))/4.0;
          if (!fuzzycmp (check))
            {
              A[vect_to_mat (m1*M2+m2,m1*M2+m2)] = 1.0;
              B[m1*M2+m2] = 0.0;
              continue;
            }

          A[vect_to_mat(m1*M2+m2,m1*M2+m2)]=check*(1.+tau/h1*fabs(V1.val(m1,m2))+tau/h2*fabs(V2.val(m1,m2)))
                                            +tau*MIU*(8./3./h1/h1+2./h2/h2);

          A[vect_to_mat(m1*M2+m2,(m1-1)*M2+m2)]=tau/2./h2*xmabs(V2.val(m1,m2))*check
                                                -tau*MIU/h2/h2;

          A[vect_to_mat(m1*M2+m2,(m1+1)*M2+m2)]=-tau/2./h2*xmabs(V2.val(m1,m2))*check
                                                -tau*MIU/h2/h2;

          A[vect_to_mat(m1*M2+m2,m1*M2+m2-1)]=tau/2./h1*xmabs(V1.val(m1,m2))*check
                                              -4./3.*tau*MIU/h1/h1;

          A[vect_to_mat(m1*M2+m2,m1*M2+m2+1)]=-tau/2./h1*xpabs(V1.val(m1,m2))*check
                                              -4./3.*tau*MIU/h1/h1;

          B[m1*M2+m2]=check*V2.val(m1,m2)-tau/h2*(p(H.right(m1,m2))-p(H.right(m1,m2-1)))
                      +tau*MIU/12./h1/h2*(V1.val(m1+1,m2+1)-V1.val(m1+1,m2-1)-V1.val(m1-1,m2+1)+V1.val(m1-1,m2-1))
                      +tau*f_3(t,x,y)*check;
        }
    }
}
}
