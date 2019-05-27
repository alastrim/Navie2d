#include "misc.h"
#include "grid.h"
#include "discrete_function.h"
#include "printers.h"
#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "loop.h"

double MIU;

int main (int argc, char **argv)
{
  std::vector<double> miu_vals = {0.1, 0.01, 0.001};
  for (double miu : miu_vals)
    {
      MIU = miu;

      std::unique_ptr<mesh> main_mesh = fillers::fill_mesh_by_arguments (argc, argv);

      timed_discrete_function tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get (), "H");
      timed_discrete_function tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V1");
      timed_discrete_function tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V2");
      trio essential (tdfH, tdfV1, tdfV2);
      fillers::fill_initial_info (essential);

      timed_discrete_function real_tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get (), "H");
      timed_discrete_function real_tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V1");
      timed_discrete_function real_tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V2");
      trio real (real_tdfH, real_tdfV1, real_tdfV2);
      fillers::fill_real_info (real);

      time_loop (essential, real);

      if (PRINT_RESULTS)
        {
          print_to_gnuplot (tdfH, "H");
          print_to_gnuplot (tdfV1, "V1");
          print_to_gnuplot (tdfV2, "V2");

          print_to_gnuplot (real_tdfH, "realH");
          print_to_gnuplot (real_tdfV1, "realV1");
          print_to_gnuplot (real_tdfV2, "realV2");
        }

      std::string filename = get_parameters_string (essential);
      print_residuals (filename, tdfH, real_tdfH, "H");
      print_residuals (filename, tdfV1, real_tdfV1, "V1");
      print_residuals (filename, tdfV2, real_tdfV2, "V2");
    }

  return 0;
}
