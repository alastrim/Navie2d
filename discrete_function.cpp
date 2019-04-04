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
  for (unsigned int i = 0; i < m_i_size; i++)
    {
      for (unsigned int j = 0; j < m_j_size; j++)
        {
          index ij = {i, j};
          point xy = m_grid->get_point (ij);
          double value = cf (xy);
          set_value (ij, value);
        }
    }
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

///////////////////////////////////////////////////////////////////////////////

timed_discrete_function::timed_discrete_function (const grid *grid, const scale *scale) : m_grid (grid), m_scale (scale)
{
  m_k_size = sbsc (m_scale->get_parameters ().m_t_step_count);

  m_data.resize (m_k_size);
}

void timed_discrete_function::fill (timed_continuous_function tcf)
{
  for (int k = 0; k < toi (m_k_size); k++)
    {
      double t = m_scale->get_time (k);
      continuous_function cf = [&] (point xy) { return tcf (t, xy); };
      std::unique_ptr<discrete_function> df = std::make_unique<discrete_function> (m_grid);
      df->fill (cf);
      set_cut (k, std::move (df));
    }
}

void timed_discrete_function::set_cut (int k, std::unique_ptr<discrete_function> df)
{
  assert (k >= 0 && tou (k) < m_k_size, "Bad argument for timed function cut");

  std::unique_ptr<discrete_function> &cut = m_data[tou (k)];
  cut = std::move (df);
}

const discrete_function &timed_discrete_function::get_cut (int k)
{
  assert (k >= 0 && tou (k) < m_k_size, "Bad argument for timed function cut");

  return *m_data[tou (k)];
}
