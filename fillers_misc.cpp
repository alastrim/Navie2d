#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "misc.h"
#include "matvec.h"
#include "discrete_function.h"
#include "grid.h"

namespace fillers
{
void fill_initial_info (trio &essential)
{
  discrete_function &H_initial_cut = essential.m_tdfH.get_cut (0);
  discrete_function &V1_initial_cut = essential.m_tdfV1.get_cut (0);
  discrete_function &V2_initial_cut = essential.m_tdfV2.get_cut (0);

  continuous_function H_initial_filler = [] (point xy) { double x = xy.first, y = xy.second; return r (0, x, y); };
  continuous_function V1_initial_filler = [] (point xy) { double x = xy.first, y = xy.second; return u1 (0, x, y); };
  continuous_function V2_initial_filler = [] (point xy) { double x = xy.first, y = xy.second; return u2 (0, x, y); };

  H_initial_cut.fill (H_initial_filler);
  V1_initial_cut.fill (V1_initial_filler);
  V2_initial_cut.fill (V2_initial_filler);
}

void fill_real_info (trio &real)
{
  real.m_tdfH.fill ([] (double t, point xy) { double x = xy.first, y = xy.second; return r (t, x, y); });
  real.m_tdfV1.fill ([] (double t, point xy) { double x = xy.first, y = xy.second; return u1 (t, x, y); });
  real.m_tdfV2.fill ([] (double t, point xy) { double x = xy.first, y = xy.second; return u2 (t, x, y); });
}

std::unique_ptr<mesh> fill_mesh_by_arguments (int argc, char **argv)
{
  int t_step_count, x_step_count, y_step_count, iT;
  if (argc < 5 || !(iT = atoi (argv[1])) || !(t_step_count = atoi(argv[2])) || !(x_step_count = atoi(argv[3])) || !(y_step_count = atoi(argv[4])))
    {
      printf ("Usage: ./main.exe <T> <t_step_count> <x_step_count> <y_step_count>\n");
      printf ("Bad arguments given, using default...\n");
      t_step_count = 35;
      x_step_count = 3;
      y_step_count = 3;
      iT = 1;
    }
  double T = iT;
  double X = 3. * M_PI;
  double Y = 3. * M_PI;

  std::unique_ptr<mesh> result = std::make_unique<mesh> ();

  grid_parameters V_grid_parameters = construct_V_grid_parameters (0, 0, X, Y, 0., 0., 2. * M_PI, 1. * M_PI, x_step_count, y_step_count);
  grid_parameters H_grid_parameters = construct_H_grid_parameters (V_grid_parameters, 0., 0., 2. * M_PI, 1. * M_PI);

  result->m_H_grid = std::make_unique<grid> (H_grid_parameters);
  result->m_V_grid = std::make_unique<grid> (V_grid_parameters);
  result->m_scale = std::make_unique<scale> (0, T, t_step_count);

  return result;
}

grid_parameters construct_H_grid_parameters (const grid_parameters &V_grid_parameters,
                                             double x_hole_origin, double y_hole_origin,
                                             double x_hole_end, double y_hole_end)
{
  grid_parameters H_grid_parameters;

  H_grid_parameters.m_x_step = V_grid_parameters.m_x_step;
  H_grid_parameters.m_y_step = V_grid_parameters.m_y_step;

  H_grid_parameters.m_x_origin = V_grid_parameters.m_x_origin + H_grid_parameters.m_x_step / 2.;
  H_grid_parameters.m_y_origin = V_grid_parameters.m_y_origin + H_grid_parameters.m_y_step / 2.;

  H_grid_parameters.m_x_point_count = V_grid_parameters.m_x_point_count - 1;
  H_grid_parameters.m_y_point_count = V_grid_parameters.m_y_point_count - 1;

  H_grid_parameters.m_hole_origin_index_x = V_grid_parameters.m_hole_origin_index_x;
  H_grid_parameters.m_hole_end_index_x = V_grid_parameters.m_hole_end_index_x;
  H_grid_parameters.m_hole_origin_index_y = V_grid_parameters.m_hole_origin_index_y;
  H_grid_parameters.m_hole_end_index_y = V_grid_parameters.m_hole_end_index_y;

  double calc_x_hole_origin = V_grid_parameters.m_x_origin + V_grid_parameters.m_hole_origin_index_x * V_grid_parameters.m_x_step;
  assert (!fuzzycmp (calc_x_hole_origin, x_hole_origin), "Hole does not conform to the grid");
  double calc_x_hole_end = V_grid_parameters.m_x_origin + V_grid_parameters.m_hole_end_index_x * V_grid_parameters.m_x_step;
  assert (!fuzzycmp (calc_x_hole_end, x_hole_end), "Hole does not conform to the grid");
  double calc_y_hole_origin = V_grid_parameters.m_y_origin + V_grid_parameters.m_hole_origin_index_y * V_grid_parameters.m_y_step;
  assert (!fuzzycmp (calc_y_hole_origin, y_hole_origin), "Hole does not conform to the grid");
  double calc_y_hole_end = V_grid_parameters.m_y_origin + V_grid_parameters.m_hole_end_index_y * V_grid_parameters.m_y_step;
  assert (!fuzzycmp (calc_y_hole_end, y_hole_end), "Hole does not conform to the grid");

  return H_grid_parameters;
}

grid_parameters construct_V_grid_parameters (double x_origin, double y_origin,
                                             double x_end, double y_end,
                                             double x_hole_origin, double y_hole_origin,
                                             double x_hole_end, double y_hole_end,
                                             int x_step_count, int y_step_count)
{
  grid_parameters V_grid_parameters;
  assert (x_step_count > 0 && y_step_count > 0, "Bad arguments for grid creation");
  assert (x_end > x_origin && y_end > y_origin, "Bad arguments for grid creation");

  V_grid_parameters.m_x_origin = x_origin;
  V_grid_parameters.m_y_origin = y_origin;
  V_grid_parameters.m_x_step = (x_end - x_origin) / x_step_count;
  V_grid_parameters.m_y_step = (y_end - y_origin) / y_step_count;
  V_grid_parameters.m_x_point_count = x_step_count + 1;
  V_grid_parameters.m_y_point_count = y_step_count + 1;

  for (int i = 0; i < V_grid_parameters.m_x_point_count; i++)
    {
      if (V_grid_parameters.m_x_origin + i * V_grid_parameters.m_x_step >= x_hole_origin && V_grid_parameters.m_hole_origin_index_x < 0)
        V_grid_parameters.m_hole_origin_index_x = i;
      if (V_grid_parameters.m_x_origin + i * V_grid_parameters.m_x_step > x_hole_end && V_grid_parameters.m_hole_end_index_x < 0)
        V_grid_parameters.m_hole_end_index_x = i - 1;
    }
  for (int j = 0; j < V_grid_parameters.m_y_point_count; j++)
    {
      if (V_grid_parameters.m_y_origin + j * V_grid_parameters.m_y_step >= y_hole_origin && V_grid_parameters.m_hole_origin_index_y < 0)
        V_grid_parameters.m_hole_origin_index_y = j;
      if (V_grid_parameters.m_y_origin + j * V_grid_parameters.m_y_step > y_hole_end && V_grid_parameters.m_hole_end_index_y < 0)
        V_grid_parameters.m_hole_end_index_y = j - 1;
    }
  assert (V_grid_parameters.m_hole_origin_index_x >= 0
          && V_grid_parameters.m_hole_end_index_x >= 0
          && V_grid_parameters.m_hole_origin_index_y >= 0
          && V_grid_parameters.m_hole_end_index_y >= 0,
          "Bad hole coordinates");

  double calc_x_hole_origin = V_grid_parameters.m_x_origin + V_grid_parameters.m_hole_origin_index_x * V_grid_parameters.m_x_step;
  assert (!fuzzycmp (calc_x_hole_origin, x_hole_origin), "Hole does not conform to the grid");
  double calc_x_hole_end = V_grid_parameters.m_x_origin + V_grid_parameters.m_hole_end_index_x * V_grid_parameters.m_x_step;
  assert (!fuzzycmp (calc_x_hole_end, x_hole_end), "Hole does not conform to the grid");
  double calc_y_hole_origin = V_grid_parameters.m_y_origin + V_grid_parameters.m_hole_origin_index_y * V_grid_parameters.m_y_step;
  assert (!fuzzycmp (calc_y_hole_origin, y_hole_origin), "Hole does not conform to the grid");
  double calc_y_hole_end = V_grid_parameters.m_y_origin + V_grid_parameters.m_hole_end_index_y * V_grid_parameters.m_y_step;
  assert (!fuzzycmp (calc_y_hole_end, y_hole_end), "Hole does not conform to the grid");

  return V_grid_parameters;
}
}
