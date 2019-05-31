#pragma once
#include "misc.h"

enum class point_type
{
  inner,
  edge,
  outer,
  INVALID
};

struct grid_parameters
{
  double m_x_origin;
  double m_y_origin;
  double m_x_step;
  double m_y_step;
  int m_x_point_count = -1;
  int m_y_point_count = -1;
  int m_hole_origin_index_x = -1;
  int m_hole_end_index_x = -1;
  int m_hole_origin_index_y = -1;
  int m_hole_end_index_y = -1;
};

class grid
{
public:
  grid (const grid_parameters &parameters);
  point get_point (index ij) const;
  point_type get_full_type (index ij) const;
  point_type get_grid_only_type (index ij) const;
  grid_parameters get_parameters () const;
private:
  grid_parameters m_parameters;
};

///////////////////////////////////////////////////////////////////////////////

struct scale_parameters
{
  double m_t_origin;
  double m_t_step;
  int m_t_point_count;
};

class scale
{
public:
  scale (double t_origin, double t_end, int t_step_count);
  double get_time (int k) const;
  scale_parameters get_parameters () const;
private:
  scale_parameters m_parameters;
};
