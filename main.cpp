#include "misc.h"
#include "grid.h"
#include "discrete_function.h"
#include "printers.h"
#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "loop.h"

int main (int argc, char **argv)
{
  std::unique_ptr<mesh> main_mesh = fillers::fill_mesh_by_arguments (argc, argv);

  timed_discrete_function tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get ());
  timed_discrete_function tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get ());
  timed_discrete_function tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get ());
  trio essential (tdfH, tdfV1, tdfV2);
  fillers::fill_initial_info (essential);

  time_loop (essential);

  timed_discrete_function real_tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get ());
  timed_discrete_function real_tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get ());
  timed_discrete_function real_tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get ());
  trio real (real_tdfH, real_tdfV1, real_tdfV2);
  fillers::fill_real_info (real);

  print_to_gnuplot (tdfH, "H");
  print_to_gnuplot (tdfV1, "V1");
  print_to_gnuplot (tdfV2, "V2");

  print_to_gnuplot (real_tdfH, "realH");
  print_to_gnuplot (real_tdfV1, "realV1");
  print_to_gnuplot (real_tdfV2, "realV2");
  return 0;
}
