#include "grid.h"

grid::grid(double x_origin, double y_origin, double x_end, double y_end, int x_step_count, int y_step_count)
{
  assert (x_step_count > 0 && y_step_count > 0, "Bad arguments for grid creation");
  assert (x_end > x_origin && y_end > y_origin, "Bad arguments for grid creation");

  m_parameters.m_x_origin = x_origin;
  m_parameters.m_y_origin = y_origin;
  m_parameters.m_x_step = (x_end - x_origin) / x_step_count;
  m_parameters.m_y_step = (y_end - y_origin) / y_step_count;
  m_parameters.m_x_step_count = x_step_count;
  m_parameters.m_y_step_count = y_step_count;
}

point grid::get_point (index ij) const
{
  assert (ij.first >= 0 && ij.first <= m_parameters.m_x_step_count, "Bad index for point");
  assert (ij.second >= 0 && ij.second <= m_parameters.m_y_step_count, "Bad index for point");

  double x = m_parameters.m_x_origin + m_parameters.m_x_step * ij.first;
  double y = m_parameters.m_y_origin + m_parameters.m_y_step * ij.second;
  return {x, y};
}

grid_parameters grid::get_parameters () const
{
  return m_parameters;
}

///////////////////////////////////////////////////////////////////////////////

scale::scale (double t_origin, double t_end, int t_step_count)
{
  assert (t_step_count > 0, "Bad arguments for scale creation");
  assert (t_end > t_origin, "Bad arguments for scale creation");

  m_parameters.m_t_origin = t_origin;
  m_parameters.m_t_step = (t_end - t_origin) / t_step_count;
  m_parameters.m_t_step_count = t_step_count;
}

double scale::get_time (int k) const
{
  assert (k >= 0 && k <= m_parameters.m_t_step_count, "Bad index for point");

  double t = m_parameters.m_t_origin + m_parameters.m_t_step * k;
  return t;
}

scale_parameters scale::get_parameters () const
{
  return m_parameters;
}
