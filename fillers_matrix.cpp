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

class Setter
{
public:
  Setter (std::vector<double> &A, int M1, int M2) : A (A), M1 (M1), M2 (M2), S (M1 * M2) {}
  double & operator () (int i, int m1, int m2)
  {
    if (m1 < 0 || m1 > M1 - 1 || m2 < 0 || m2 > M2 - 1)
      return dummy;
    return A[tou (i * S + m1 * M2 + m2)];
  }
private:
  std::vector<double> &A;
  int M1;
  int M2;
  int S;
  double dummy;
};

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

  int M1 = g_parameters.m_x_step_count + 1;
  int M2 = g_parameters.m_y_step_count + 1;
  int S = static_cast<int> (H.get_raw_vector ().size ());
  assert (S == M1 * M2, "Bad matrix size");
  A = std::vector<double> (tou (S * S), 0);
  B = std::vector<double> (tou (S), 0);

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  Setter A_at (A, M1, M2);

  for (int m1 = 0; m1 < M1; m1++)
    {
      for (int m2 = 0; m2 < M2; m2++)
        {
          int i = m1 * M2 + m2;
          double t = sc->get_time (k + 1), x = gr->get_point ({m1,m2}).first, y = gr->get_point ({m1,m2}).second;

          A_at(i,m1,m2)=1.
                        +tau*(xpabs(V1.tilda(m1+1,m2))-xmabs(V1.tilda(m1,m2)))/2.0/h1
                        +tau*(xpabs(V2.tilda(m1,m2+1))-xmabs(V2.tilda(m1,m2)))/2.0/h2;

          A_at(i,m1,m2+1)=tau*xmabs(V2.tilda (m1,m2+1))/2.0/h2;
          A_at(i,m1+1,m2)=tau*xmabs(V1.tilda (m1+1,m2))/2.0/h1;
          A_at(i,m1-1,m2)=-tau*xpabs(V1.tilda (m1,m2))/2.0/h1;
          A_at(i,m1,m2-1)=-tau*xpabs(V2.tilda (m1,m2))/2.0/h2;

          B[tou(i)]=H.val(m1,m2)+tau*f_1(t,x,y);
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

  int M1 = g_parameters.m_x_step_count + 1;
  int M2 = g_parameters.m_y_step_count + 1;
  int S = static_cast<int> (V1.get_raw_vector ().size ());
  assert (S == (M1) * (M2), "Bad matrix size");
  A = std::vector<double> (tou(S * S), 0);
  B = std::vector<double> (tou(S), 0);

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  Setter A_at (A, M1, M2);

  for (int m1 = 0; m1 < M1; m1++)
    {
      for (int m2 = 0; m2 < M2; m2++)
        {
          int i = m1 * M2 + m2;
          double t = sc->get_time (k + 1), x = gr->get_point ({m1,m2}).first, y = gr->get_point ({m1,m2}).second;
          if (m1 == 0 || m2 == 0 || m1 == (M1 - 1) || m2 == (M2 - 1))
            {
              A_at(i,m1,m2) = 1.0;
              B[tou(i)] = 0;
              continue;
            }

          double check = (H.val(m1,m2)+H.val(m1,m2-1)+H.val(m1-1,m2)+H.val(m1-1,m2-1))/4.0;
          if (!fuzzycmp (check))
            {
              A_at(i,m1,m2) = 1.0;
              B[tou(i)] = 0;
              continue;
            }

          A_at(i,m1,m2)=check*(1.+tau/h1*fabs(V1.val(m1,m2))+tau/h2*fabs(V2.val(m1,m2)))
                        +tau*MIU*(8./3./h1/h1+2./h2/h2);

          A_at(i,m1-1,m2)=-tau/2./h1*xpabs(V1.val(m1,m2))*check
                          -4./3.*tau*MIU/h1/h1;

          A_at(i,m1+1,m2)=tau/2./h1*xmabs(V1.val(m1,m2))*check
                          -4./3.*tau*MIU/h1/h1;

          A_at(i,m1,m2-1)=-tau/2./h2*xpabs(V2.val(m1,m2))*check
                          -tau*MIU/h2/h2;

          A_at(i,m1,m2+1)=tau/2./h2*xmabs(V2.val(m1,m2))*check
                          -tau*MIU/h2/h2;

          B[tou(i)]=check*V1.val(m1,m2)-tau/h1*(p(H.left(m1,m2))-p(H.left(m1-1,m2)))
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
  const grid *gr = essential.m_tdfV2.get_grid ();
  const scale *sc = essential.m_tdfV2.get_scale ();
  grid_parameters g_parameters = gr->get_parameters ();
  scale_parameters s_parameters = sc->get_parameters ();

  int M1 = g_parameters.m_x_step_count + 1;
  int M2 = g_parameters.m_y_step_count + 1;
  int S = static_cast<int> (V2.get_raw_vector ().size ());
  assert (S == (M1) * (M2), "Bad matrix size");
  A = std::vector<double> (tou(S * S), 0);
  B = std::vector<double> (tou(S), 0);

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  Setter A_at (A, M1, M2);

  for (int m1 = 0; m1 < M1; m1++)
    {
      for (int m2 = 0; m2 < M2; m2++)
        {
          int i = m1 * M2 + m2;
          double t = sc->get_time (k + 1), x = gr->get_point ({m1,m2}).first, y = gr->get_point ({m1,m2}).second;

          if (m1 == 0 || m2 == 0 || m1 == (M1 - 1) || m2 == (M2 - 1))
            {
              A_at(i,m1,m2) = 1.0;
              B[tou(i)] = 0.0;
              continue;
            }

          double check = (H.val(m1,m2)+H.val(m1,m2-1)+H.val(m1-1,m2)+H.val(m1-1,m2-1))/4.0;
          if (!fuzzycmp (check))
            {
              A_at(i,m1,m2) = 1.0;
              B[tou(i)] = 0.0;
              continue;
            }

          A_at(i,m1,m2)=check*(1.+tau/h1*fabs(V1.val(m1,m2))+tau/h2*fabs(V2.val(m1,m2)))
                        +tau*MIU*(2./h1/h1+8./3./h2/h2);

          A_at(i,m1-1,m2)=-tau/2./h1*xpabs(V1.val(m1,m2))*check
                          -tau*MIU/h1/h1;

          A_at(i,m1+1,m2)=tau/2./h1*xmabs(V1.val(m1,m2))*check
                          -tau*MIU/h1/h1;

          A_at(i,m1,m2-1)=-tau/2./h2*xpabs(V2.val(m1,m2))*check
                          -4./3.*tau*MIU/h2/h2;

          A_at(i,m1,m2+1)=tau/2./h2*xmabs(V2.val(m1,m2))*check
                          -4./3.*tau*MIU/h2/h2;

          B[tou(i)]=check*V2.val(m1,m2)-tau/h2*(p(H.right(m1,m2))-p(H.right(m1,m2-1)))
                    +tau*MIU/12./h1/h2*(V1.val(m1+1,m2+1)-V1.val(m1+1,m2-1)-V1.val(m1-1,m2+1)+V1.val(m1-1,m2-1))
                    +tau*f_3(t,x,y)*check;
        }
    }
}
}
