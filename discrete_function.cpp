#include "discrete_function.h"
#include "grid.h"

discrete_function::discrete_function (const grid *grid) : m_grid (grid)
{
  m_i_size = sbsc (m_grid->get_parameters ().m_x_step_count);
  m_j_size = sbsc (m_grid->get_parameters ().m_y_step_count);

  m_data.resize (m_i_size * m_j_size);
}

void discrete_function::fill (continuous_function cf)
{
  do_for_each ([&cf] (index ij, point xy, discrete_function &self){
    double value = cf (xy);
    self.set_value (ij, value);
  });
}

void discrete_function::set_value (index ij, double value)
{
  assert (ij.first >= 0 && tou (ij.first) < m_i_size, "Bad argument for function set value");
  assert (ij.second >= 0 && tou (ij.second) < m_j_size, "Bad argument for function set value");

  m_data[tou (ij.first) * m_j_size + tou (ij.second)] = value;
}

double discrete_function::get_value (index ij)
{
  assert (ij.first >= 0 && tou (ij.first) < m_i_size, "Bad argument for function get value");
  assert (ij.second >= 0 && tou (ij.second) < m_j_size, "Bad argument for function get value");

  return m_data[tou (ij.first) * m_j_size + tou (ij.second)];
}

void discrete_function::do_for_each (discrete_foreach_function dff)
{
  for (unsigned int i = 0; i < m_i_size; i++)
    {
      for (unsigned int j = 0; j < m_j_size; j++)
        {
          index ij = {i, j};
          point xy = m_grid->get_point (ij);
          dff (ij, xy, *this);
        }
    }
}

void discrete_function::do_for_edge (discrete_foreach_function dff)
{
  index ij;
  point xy;

  for (unsigned int i = 0; i < m_i_size; i++)
    {
      ij = {i, 0};
      xy = m_grid->get_point (ij);
      dff (ij, xy, *this);

      ij = {i, m_j_size - 1};
      xy = m_grid->get_point (ij);
      dff (ij, xy, *this);
    }

  for (unsigned int j = 1; j < m_j_size - 1; j++)
    {
      ij = {0, j};
      xy = m_grid->get_point (ij);
      dff (ij, xy, *this);

      ij = {m_i_size - 1, j};
      xy = m_grid->get_point (ij);
      dff (ij, xy, *this);
    }
}

///////////////////////////////////////////////////////////////////////////////

timed_discrete_function::timed_discrete_function (const grid *grid, const scale *scale) : m_grid (grid), m_scale (scale)
{
  m_k_size = sbsc (m_scale->get_parameters ().m_t_step_count);

  m_data.resize (m_k_size);
}

void timed_discrete_function::fill (timed_continuous_function tcf)
{
  do_for_each ([&tcf] (int k, double t, timed_discrete_function &self) {
      continuous_function cf = [&] (point xy) { return tcf (t, xy); };
      std::unique_ptr<discrete_function> df = std::make_unique<discrete_function> (self.m_grid);
      df->fill (cf);
      self.set_cut (k, std::move (df));
    });
}

void timed_discrete_function::set_cut (int k, std::unique_ptr<discrete_function> df)
{
  assert (k >= 0 && tou (k) < m_k_size, "Bad argument for timed function cut");

  std::unique_ptr<discrete_function> &cut = m_data[tou (k)];
  cut = std::move (df);
}

discrete_function &timed_discrete_function::get_cut (int k)
{
  assert (k >= 0 && tou (k) < m_k_size, "Bad argument for timed function cut");

  return *m_data[tou (k)];
}

void timed_discrete_function::do_for_each (timed_discrete_foreach_function tdff)
{
  for (int k = 0; k < toi (m_k_size); k++)
    {
      double t = m_scale->get_time (k);
      tdff (k, t, *this);
    }
}

