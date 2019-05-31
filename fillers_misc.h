#pragma once
#include "misc.h"
struct grid_parameters;

struct mesh
{
  std::unique_ptr<grid> m_H_grid;
  std::unique_ptr<grid> m_V_grid;
  std::unique_ptr<scale> m_scale;
};

struct trio
{
  trio (const trio &) = delete;
  trio (timed_discrete_function &tdfH, timed_discrete_function &tdfV1, timed_discrete_function &tdfV2)
    : m_tdfH (tdfH), m_tdfV1 (tdfV1), m_tdfV2 (tdfV2) {}
  timed_discrete_function &m_tdfH;
  timed_discrete_function &m_tdfV1;
  timed_discrete_function &m_tdfV2;
};

namespace fillers
{
grid_parameters construct_V_grid_parameters (double x_origin, double y_origin,
                                             double x_end, double y_end,
                                             double x_hole_origin, double y_hole_origin,
                                             double x_hole_end, double y_hole_end,
                                             int x_step_count, int y_step_count);
grid_parameters construct_H_grid_parameters (const grid_parameters &V_grid_parameters,
                                             double x_hole_origin, double y_hole_origin,
                                             double x_hole_end, double y_hole_end);
void fill_initial_info (trio &essential);
void fill_real_info (trio &real);
std::unique_ptr<mesh> fill_mesh_by_arguments (int argc, char **argv);
}
