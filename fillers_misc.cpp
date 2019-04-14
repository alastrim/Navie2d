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

  discrete_foreach_function zero_setter = [&] (index ij, point, discrete_function &self) { self.set_value (ij, 0); };
  timed_discrete_foreach_function timed_zero_setter = [&] (int k, double, timed_discrete_function &self) { self.get_cut (k).do_for_edge (zero_setter); };

//  essential.m_tdfH.do_for_each (timed_zero_setter);
  essential.m_tdfV1.do_for_each (timed_zero_setter);
  essential.m_tdfV2.do_for_each (timed_zero_setter);
}

void fill_real_info (trio &real)
{
  real.m_tdfH.fill ([] (double t, point xy) { double x = xy.first, y = xy.second; return r (t, x, y); });
  real.m_tdfV1.fill ([] (double t, point xy) { double x = xy.first, y = xy.second; return u1 (t, x, y); });
  real.m_tdfV2.fill ([] (double t, point xy) { double x = xy.first, y = xy.second; return u2 (t, x, y); });

  discrete_foreach_function zero_setter = [&] (index ij, point, discrete_function &self) { self.set_value (ij, 0); };
  timed_discrete_foreach_function timed_zero_setter = [&] (int k, double, timed_discrete_function &self) { self.get_cut (k).do_for_edge (zero_setter); };

//  real.m_tdfH.do_for_each (timed_zero_setter);
  real.m_tdfV1.do_for_each (timed_zero_setter);
  real.m_tdfV2.do_for_each (timed_zero_setter);
}

std::unique_ptr<mesh> fill_mesh_by_arguments (int argc, char **argv)
{
  int t_step_count, x_step_count, y_step_count, iT, iX, iY;
  if (argc < 7 || !(iT = atoi (argv[1])) || !(t_step_count = atoi(argv[2])) || !(iX = atoi (argv[3])) || !(x_step_count = atoi(argv[4])) || !(iY = atoi (argv[5])) || !(y_step_count = atoi(argv[6])))
    {
      printf ("Usage: ./main.exe <T> <t_step_count> <X> <x_step_count> <Y> <y_step_count>\n");
      printf ("Bad arguments given, using default...\n");
      t_step_count = 9;
      x_step_count = 3;
      y_step_count = 3;
      iT = 1;
      iX = 1;
      iY = 1;
    }
  double T = iT, X = iX, Y = iY;

  std::unique_ptr<mesh> result = std::make_unique<mesh> ();

  double half_x_step = (X - 0) / x_step_count / 2;
  double half_y_step = (Y - 0) / y_step_count / 2;

  result->m_H_grid = std::make_unique<grid> (0 + half_x_step, 0 + half_y_step,
                                             X - half_x_step, Y - half_y_step,
                                             x_step_count - 1, y_step_count - 1);
  result->m_V_grid = std::make_unique<grid> (0, 0, X, Y, x_step_count, y_step_count);
  result->m_scale = std::make_unique<scale> (0, T, t_step_count);

  return result;
}
}
