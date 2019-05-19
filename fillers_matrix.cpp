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


class MatrixSetter
{
public:
  MatrixSetter (std::vector<double> &A, const discrete_function &df) : A (A), m_df (df) {}
  double & operator () (int m1_base, int m2_base, int m1_mod, int m2_mod)
  {
    const grid *gr = m_df.get_grid ();
    grid_parameters gr_p = gr->get_parameters ();
    int M1 = gr_p.m_x_point_count;
    int M2 = gr_p.m_y_point_count;
    int S = M1 * M2;

    int i = m1_base * M2 + m2_base;
    int m1 = m1_base + m1_mod, m2 = m2_base + m2_mod;


    if (gr->get_type ({m1, m2}) == point_type::outer)
      return m_dummy;

    return A[tou (i * S + m1 * M2 + m2)];
  }
private:
  std::vector<double> &A;
  const discrete_function &m_df;
  double m_dummy;
};

class VectorSetter
{
public:
  VectorSetter (std::vector<double> &B, const discrete_function &df) : B (B), m_df (df) {}
  double & operator () (int m1, int m2)
  {
    const grid *gr = m_df.get_grid ();
    grid_parameters gr_p = gr->get_parameters ();
    int M2 = gr_p.m_y_point_count;

    return B[tou (m1 * M2 + m2)];
  }
private:
  std::vector<double> &B;
  const discrete_function &m_df;
};

static void trivial (int m1, int m2, MatrixSetter &A_at, VectorSetter &B_at)
{
  A_at (m1,m2,0,0) = 1.0;
  B_at (m1,m2) = 0.0;
}

namespace fillers
{
void fill_first (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
  discrete_function &H = essential.m_tdfH.get_cut (k), &V1 = essential.m_tdfV1.get_cut (k), &V2 = essential.m_tdfV2.get_cut (k);
  const grid *gr = essential.m_tdfH.get_grid ();
  const scale *sc = essential.m_tdfH.get_scale ();
  grid_parameters g_parameters = gr->get_parameters ();
  scale_parameters s_parameters = sc->get_parameters ();

  double tau = s_parameters.m_t_step, h1 = g_parameters.m_x_step, h2 = g_parameters.m_y_step;

  MatrixSetter A_at (A, H);
  VectorSetter B_at (B, H);

  H.do_for_each ([&] (index ij, point xy, discrete_function &)
  {
    int m1 = ij.first, m2 = ij.second;
    double t = sc->get_time (k + 1), x = xy.first, y = xy.second;

    A_at(m1,m2,0,0)=1.
                  +tau*(xpabs(V1.tilda(m1+1,m2))-xmabs(V1.tilda(m1,m2)))/2.0/h1
                  +tau*(xpabs(V2.tilda(m1,m2+1))-xmabs(V2.tilda(m1,m2)))/2.0/h2;

    A_at(m1,m2,0,+1)=tau*xmabs(V2.tilda (m1,m2+1))/2.0/h2;
    A_at(m1,m2,+1,0)=tau*xmabs(V1.tilda (m1+1,m2))/2.0/h1;
    A_at(m1,m2,-1,0)=-tau*xpabs(V1.tilda (m1,m2))/2.0/h1;
    A_at(m1,m2,0,-1)=-tau*xpabs(V2.tilda (m1,m2))/2.0/h2;

    B_at(m1,m2)=H.val(m1,m2)+tau*f_1(t,x,y);
  });
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

  int M1 = g_parameters.m_x_point_count;
  int M2 = g_parameters.m_y_point_count;
  int S = static_cast<int> (V1.get_raw_vector ().size ());
  assert (S == (M1) * (M2), "Bad matrix size");
  A = std::vector<double> (tou(S * S), 0);
  B = std::vector<double> (tou(S), 0);

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  MatrixSetter A_at (A, V1);
  VectorSetter B_at (B, V1);

  V1.do_for_each ([&] (index ij, point xy, discrete_function &)
  {
    int m1 = ij.first, m2 = ij.second;
    double t = sc->get_time (k + 1), x = xy.first, y = xy.second;
    double check = (H.val(m1,m2)+H.val(m1,m2-1)+H.val(m1-1,m2)+H.val(m1-1,m2-1))/4.0;

    if (gr->get_type ({m1, m2}) == point_type::edge || !fuzzycmp (check))
      return trivial (m1, m2, A_at, B_at);

    A_at(m1,m2,0,0)=check*(1.+tau/h1*fabs(V1.val(m1,m2))+tau/h2*fabs(V2.val(m1,m2)))
                    +tau*MIU*(8./3./h1/h1+2./h2/h2);

    A_at(m1,m2,-1,0)=-tau/2./h1*xpabs(V1.val(m1,m2))*check
                    -4./3.*tau*MIU/h1/h1;

    A_at(m1,m2,+1,0)=tau/2./h1*xmabs(V1.val(m1,m2))*check
                    -4./3.*tau*MIU/h1/h1;

    A_at(m1,m2,0,-1)=-tau/2./h2*xpabs(V2.val(m1,m2))*check
                    -tau*MIU/h2/h2;

    A_at(m1,m2,0,+1)=tau/2./h2*xmabs(V2.val(m1,m2))*check
                    -tau*MIU/h2/h2;

    B_at(m1,m2)=check*V1.val(m1,m2)-tau/h1*(p(H.left(m1,m2))-p(H.left(m1-1,m2)))
              +tau*MIU/12./h1/h2*(V2.val(m1+1,m2+1)-V2.val(m1+1,m2-1)-V2.val(m1-1,m2+1)+V2.val(m1-1,m2-1))
              +tau*f_2(t,x,y)*check;
  });
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

  int M1 = g_parameters.m_x_point_count;
  int M2 = g_parameters.m_y_point_count;
  int S = static_cast<int> (V2.get_raw_vector ().size ());
  assert (S == (M1) * (M2), "Bad matrix size");
  A = std::vector<double> (tou(S * S), 0);
  B = std::vector<double> (tou(S), 0);

  double tau = s_parameters.m_t_step;
  double h1 = g_parameters.m_x_step;
  double h2 = g_parameters.m_y_step;

  MatrixSetter A_at (A, V2);
  VectorSetter B_at (B, V2);

  V1.do_for_each ([&] (index ij, point xy, discrete_function &)
  {
    int m1 = ij.first, m2 = ij.second;
    double t = sc->get_time (k + 1), x = xy.first, y = xy.second;
    double check = (H.val(m1,m2)+H.val(m1,m2-1)+H.val(m1-1,m2)+H.val(m1-1,m2-1))/4.0;

    if (gr->get_type ({m1, m2}) == point_type::edge || !fuzzycmp (check))
      return trivial (m1, m2, A_at, B_at);

    A_at(m1,m2,0,0)=check*(1.+tau/h1*fabs(V1.val(m1,m2))+tau/h2*fabs(V2.val(m1,m2)))
                    +tau*MIU*(2./h1/h1+8./3./h2/h2);

    A_at(m1,m2,-1,0)=-tau/2./h1*xpabs(V1.val(m1,m2))*check
                     -tau*MIU/h1/h1;

    A_at(m1,m2,+1,0)=tau/2./h1*xmabs(V1.val(m1,m2))*check
                     -tau*MIU/h1/h1;

    A_at(m1,m2,0,-1)=-tau/2./h2*xpabs(V2.val(m1,m2))*check
                     -4./3.*tau*MIU/h2/h2;

    A_at(m1,m2,0,+1)=tau/2./h2*xmabs(V2.val(m1,m2))*check
                    -4./3.*tau*MIU/h2/h2;

    B_at(m1,m2)=check*V2.val(m1,m2)-tau/h2*(p(H.right(m1,m2))-p(H.right(m1,m2-1)))
              +tau*MIU/12./h1/h2*(V1.val(m1+1,m2+1)-V1.val(m1+1,m2-1)-V1.val(m1-1,m2+1)+V1.val(m1-1,m2-1))
              +tau*f_3(t,x,y)*check;
  });
}
}
