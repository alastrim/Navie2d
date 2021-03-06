#include "setter_helpers.h"
#include "discrete_function.h"
#include "grid.h"
#include "fillers_functions.h"

double xmabs (double x)
{
  return x - fabs (x);
}

double xpabs (double x)
{
  return x + fabs (x);
}

MatrixSetter::MatrixSetter (std::map<unsigned int, double> &A, const discrete_function &df)
  : m_A (A), m_df (df)
{}

VectorSetter::VectorSetter (std::vector<double> &B, const discrete_function &df)
  : m_B (B), m_df (df)
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


  if (gr->get_grid_only_type ({m1, m2}) == point_type::outer)
    return m_dummy;

  int ind = i * S + m1 * M2 + m2;
  assert (std::find (m_taken.begin (), m_taken.end (), ind) == m_taken.end (), "Taken sanity");
  m_taken.push_back (ind);

  return m_A[tou (ind)];
}

double & VectorSetter::operator () (int m1, int m2)
{
  const grid *gr = m_df.get_grid ();
  grid_parameters gr_p = gr->get_parameters ();
  int M2 = gr_p.m_y_point_count;

  int ind = m1 * M2 + m2;
  assert (std::find (m_taken.begin (), m_taken.end (), ind) == m_taken.end (), "Taken sanity");
  m_taken.push_back (ind);

  return m_B[tou (ind)];
}

bool process_V_condition (int m1, int m2, MatrixSetter &A_at, VectorSetter &B_at)
{
  if (KNOWN_FUNC)
    return false;
  const discrete_function &V2 = A_at.m_df;
  if (m2 == 0)
    {
      B_at (m1, m2) = V2.val (m1, m2 + 1);
      A_at (m1, m2, 0, 0) = 1.;
      return true;
    }
  return false;
}

bool process_V_edge (int m1, int m2, double check, MatrixSetter &A_at, VectorSetter &B_at)
{
  const grid *gr = A_at.m_df.get_grid ();
  if (gr->get_full_type ({m1, m2}) == point_type::outer || gr->get_full_type ({m1, m2}) == point_type::edge || !fuzzycmp (check))
    {
      A_at (m1,m2,0,0) = 1.0;
      B_at (m1,m2) = 0.0;
      return true;
    }
  return false;
}

bool process_H_edge (int m1, int m2, MatrixSetter &A_at, VectorSetter &B_at)
{
  const grid *gr = A_at.m_df.get_grid ();
  if (gr->get_full_type ({m1, m2}) == point_type::outer)
    {
      A_at (m1,m2,0,0) = 1.0;
      B_at (m1,m2) = 0.0;
      return true;
    }
  return false;
}
