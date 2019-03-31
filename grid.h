#ifndef GRID_H
#define GRID_H

#include "misc.h"

class grid
{
public:
  grid();
private:
  double m_x_origin;
  double m_y_origin;
  double m_x_step;
  double m_y_step;
  int m_x_step_count;
  int m_y_step_count;
};

#endif // GRID_H
