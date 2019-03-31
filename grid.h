#ifndef GRID_H
#define GRID_H

#include "misc.h"

typedef std::pair<double, double> point;
typedef std::pair<int, int> index;

struct grid_parameters
{
  double m_x_origin;
  double m_y_origin;
  double m_x_step;
  double m_y_step;
  int m_x_step_count;
  int m_y_step_count;
};

class grid
{
public:
  grid (double x_origin, double y_origin, double x_end, double y_end, int x_step_count, int y_step_count);
  point get_point (index ij);
  grid_parameters get_parameters ();
private:
  grid_parameters m_parameters;
};

///////////////////////////////////////////////////////////////////////////////

struct scale_parameters
{
  double m_t_origin;
  double m_t_step;
  int m_t_step_count;
};

class scale
{
public:
  scale (double t_origin, double t_end, int t_step_count);
  double get_time (int k);
  scale_parameters get_parameters ();
private:
  scale_parameters m_parameters;
};

#endif // GRID_H
