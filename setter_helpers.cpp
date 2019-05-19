#include "setter_helpers.h"
#include "discrete_function.h"
#include "grid.h"

double xmabs (double x)
{
  return x - fabs (x);
}

double xpabs (double x)
{
  return x + fabs (x);
}

MatrixSetter::MatrixSetter (std::vector<double> &A, const discrete_function &df)
  : A (A), m_df (df)
{}

VectorSetter::VectorSetter (std::vector<double> &B, const discrete_function &df)
  : B (B), m_df (df)
{}

double & MatrixSetter::operator () (int m1_base, int m2_base, int m1_mod, int m2_mod)
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

double & VectorSetter::operator () (int m1, int m2)
{
  const grid *gr = m_df.get_grid ();
  grid_parameters gr_p = gr->get_parameters ();
  int M2 = gr_p.m_y_point_count;

  return B[tou (m1 * M2 + m2)];
}

bool process_if_edge (int m1, int m2, double check, MatrixSetter &A_at, VectorSetter &B_at)
{
  const grid *gr = A_at.m_df.get_grid ();
  if (gr->get_type ({m1, m2}) == point_type::edge || !fuzzycmp (check))
    {
      A_at (m1,m2,0,0) = 1.0;
      B_at (m1,m2) = 0.0;
      return true;
    }
  return false;
}
