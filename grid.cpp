#include "grid.h"

grid::grid (const grid_parameters &parameters) : m_parameters (parameters)
{
}

point_type grid::get_full_type (index ij) const
{
  int i = ij.first;
  int j = ij.second;

  point_type res_point_type = get_grid_only_type (ij);

  if (i == m_parameters.m_hole_end_index_x
      && j >= m_parameters.m_hole_origin_index_y && j <= m_parameters.m_hole_end_index_y)
    return point_type::edge;
  if (i >= m_parameters.m_hole_origin_index_x && i <= m_parameters.m_hole_end_index_x
      && j == m_parameters.m_hole_end_index_y)
    return point_type::edge;
  if (i >= m_parameters.m_hole_origin_index_x && i <= m_parameters.m_hole_end_index_x
      && j >= m_parameters.m_hole_origin_index_y && j <= m_parameters.m_hole_end_index_y)
    return point_type::outer;

  return res_point_type;
}

point_type grid::get_grid_only_type (index ij) const
{
  int i = ij.first;
  int j = ij.second;

  point_type x_point_type = point_type::INVALID;
  point_type y_point_type = point_type::INVALID;
  point_type res_point_type = point_type::INVALID;

  if (i == 0 || i == m_parameters.m_x_point_count - 1)
    x_point_type = point_type::edge;
  else if (i > 0 && i < m_parameters.m_x_point_count - 1)
    x_point_type = point_type::inner;
  else if (i < 0 || i > m_parameters.m_x_point_count - 1)
    x_point_type = point_type::outer;
  assert (x_point_type != point_type::INVALID, "Bad point type");

  if (j == 0 || j == m_parameters.m_y_point_count - 1)
    y_point_type = point_type::edge;
  else if (j > 0 && j < m_parameters.m_y_point_count - 1)
    y_point_type = point_type::inner;
  else if (j < 0 || j > m_parameters.m_y_point_count - 1)
    y_point_type = point_type::outer;
  assert (y_point_type != point_type::INVALID, "Bad point type");

  if (x_point_type == point_type::outer || y_point_type == point_type::outer)
    res_point_type = point_type::outer;
  else if (x_point_type == point_type::inner && y_point_type == point_type::inner)
    res_point_type = point_type::inner;
  else if (x_point_type == point_type::edge || y_point_type == point_type::edge)
    res_point_type = point_type::edge;
  assert (res_point_type != point_type::INVALID, "Bad point type");

  return res_point_type;
}

point grid::get_point (index ij) const
{
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
  m_parameters.m_t_point_count = t_step_count + 1;
}

double scale::get_time (int k) const
{
  assert (k >= 0 && k < m_parameters.m_t_point_count, "Bad index for point");

  double t = m_parameters.m_t_origin + m_parameters.m_t_step * k;
  return t;
}

scale_parameters scale::get_parameters () const
{
  return m_parameters;
}
