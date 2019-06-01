#include "misc.h"
#include "grid.h"
#include "discrete_function.h"
#include "printers.h"
#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
#include "loop.h"

double MIU;
double OMEGA;
double RHO_GAMMA;
double RHO_NULL;

enum class residual_t
{
  C,
  L2
};

static std::string ets (residual_t type)
{
  switch (type)
    {
  case residual_t::C:
    return "C   ";
  case residual_t::L2:
    return "L2  ";
    }
  return "Err";
}

#define TSIZE 4
#define TSTART 21
#define TMOD 2

int main (int argc, char **argv)
{
  OMEGA = 1;
  RHO_NULL = 0;
  RHO_GAMMA = 1;

  if (!PRINT_RESULTS)
    {
      FILE *f = fopen ("tables.txt", "w");
      fprintf (f, "STARTED ANOTHER RUN\n");
      fclose (f);
    }

  std::vector<double> miu_vals = {0.1, 0.01, 0.001};
  for (double miu_val : miu_vals)
    {
      MIU = miu_val;
      if (!PRINT_RESULTS)
        {
          assert (KNOWN_FUNC, "Why create table if you cant calculate residuals?");

          typedef std::tuple<double, double, double> table_value_t;
          typedef std::vector<std::pair<int, table_value_t>> table_line_t; // M value
          typedef std::vector<std::pair<int, table_line_t>> table_t; // N value
          typedef std::vector<std::pair<residual_t, table_t>> table_set_t; // type of residual

          std::vector<residual_t> types = {residual_t::C, residual_t::L2};
          table_set_t set;
          for (auto t : types)
            {
              table_value_t e (-1, -1, -1);
              table_line_t em (TSIZE, {0, e});
              table_t emp (TSIZE, {0, em});
              set.push_back ({t, emp});
            }

          for (unsigned int Nloop = 0; Nloop < TSIZE; Nloop++)
            {
              for (unsigned int Mloop = 0; Mloop < TSIZE; Mloop++)
                {
                  printf ("\nStarting another calculation...\n");
                  std::unique_ptr<mesh> main_mesh;
                  int N = TSTART * static_cast<int> (pow (TMOD, static_cast<double> (Nloop)));
                  int M = TSTART * static_cast<int> (pow (TMOD, static_cast<double> (Mloop)));

                  int t_step_count = N, x_step_count = M, y_step_count = M;
                  double T = 1, X = 3. * M_PI, Y = 3. * M_PI;
                  std::unique_ptr<mesh> result = std::make_unique<mesh> ();

                  grid_parameters V_grid_parameters = fillers::construct_V_grid_parameters (0, 0, X, Y, 0., 0., 2. * M_PI, 1. * M_PI, x_step_count, y_step_count);
                  grid_parameters H_grid_parameters = fillers::construct_H_grid_parameters (V_grid_parameters, 0., 0., 2. * M_PI, 1. * M_PI);

                  result->m_H_grid = std::make_unique<grid> (H_grid_parameters);
                  result->m_V_grid = std::make_unique<grid> (V_grid_parameters);
                  result->m_scale = std::make_unique<scale> (0, T, t_step_count);
                  main_mesh = std::move (result);


                  timed_discrete_function tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get (), "H");
                  timed_discrete_function tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V1");
                  timed_discrete_function tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V2");
                  trio essential (tdfH, tdfV1, tdfV2);
                  fillers::fill_initial_info (essential);
                  print_parameters (essential);

                  timed_discrete_function real_tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get (), "H");
                  timed_discrete_function real_tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V1");
                  timed_discrete_function real_tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V2");
                  trio real (real_tdfH, real_tdfV1, real_tdfV2);
                  fillers::fill_real_info (real);

                  time_loop (essential, &real);

                  double t1 = essential.m_tdfH.residual_C (real.m_tdfH);
                  double t2 = essential.m_tdfV1.residual_C (real.m_tdfV1);
                  double t3 = essential.m_tdfV2.residual_C (real.m_tdfV2);
                  double t4 = essential.m_tdfH.residual_L2 (real.m_tdfH);
                  double t5 = essential.m_tdfV1.residual_L2 (real.m_tdfV1);
                  double t6 = essential.m_tdfV2.residual_L2 (real.m_tdfV2);

                  set[0].second[Nloop].first = N;
                  set[0].second[Nloop].second[Mloop].first = M;
                  set[0].second[Nloop].second[Mloop].second = {t1, t2, t3};
                  set[1].second[Nloop].first = N;
                  set[1].second[Nloop].second[Mloop].first = M;
                  set[1].second[Nloop].second[Mloop].second = {t4, t5, t6};

                  print_residuals (tdfH, real_tdfH, "H");
                  print_residuals (tdfV1, real_tdfV1, "V1");
                  print_residuals (tdfV2, real_tdfV2, "V2");

                }
            }

          FILE *f = fopen ("tables.txt", "a");
          fprintf (f, "\n\n\n\\subsection{MIU = %.3f}\n", MIU);
          for (std::pair<residual_t, table_t> &table : set)
            {
              fprintf (f, "\\begin{table}[H]\n"
                      "\\caption {Функция: %s Тип невязки: %s}\n"
                      "\\begin{center}\n"
                      "\\begin{tabular}{l|l|l|l|l}\n"
                      "\\hline\n"
                      "M/N  & 21 & 42 & 84 & 168 \\\\ \\hline\n", "H", ets (table.first).c_str ());
              for (std::pair<int, table_line_t> &line : table.second)
                {
                  fprintf (f, "%4d ", line.first);
                  for (std::pair<int, table_value_t> & value : line.second)
                    fprintf (f, "& %1.0e", std::get<0>(value.second));
                  fprintf (f, "\\\\ \\hline\n");
                }
              fprintf (f, "\\end{tabular}\n"
                      "\\end{center}\n"
                      "\\end{table}\n");

              fprintf (f, "\\begin{table}[H]\n"
                          "\\caption {Функция: %s Тип невязки: %s}\n"
                      "\\begin{center}\n"
                      "\\begin{tabular}{l|l|l|l|l}\n"
                      "\\hline\n"
                      "M/N  & 21 & 42 & 84 & 168 \\\\ \\hline\n", "V1", ets (table.first).c_str ());
              for (std::pair<int, table_line_t> &line : table.second)
                {
                  fprintf (f, "%4d ", line.first);
                  for (std::pair<int, table_value_t> & value : line.second)
                    fprintf (f, "& %1.0e", std::get<1>(value.second));
                  fprintf (f, "\\\\ \\hline\n");
                }
              fprintf (f, "\\end{tabular}\n"
                      "\\end{center}\n"
                      "\\end{table}\n");

              fprintf (f, "\\begin{table}[H]\n"
                          "\\caption {Функция: %s Тип невязки: %s}\n"
                      "\\begin{center}\n"
                      "\\begin{tabular}{l|l|l|l|l}\n"
                      "\\hline\n"
                      "M/N  & 21 & 42 & 84 & 168 \\\\ \\hline\n", "V2", ets (table.first).c_str ());
              for (std::pair<int, table_line_t> &line : table.second)
                {
                  fprintf (f, "%4d ", line.first);
                  for (std::pair<int, table_value_t> & value : line.second)
                    fprintf (f, "& %1.0e", std::get<2>(value.second));
                  fprintf (f, "\\\\ \\hline\n");
                }
              fprintf (f, "\\end{tabular}\n"
                      "\\end{center}\n"
                      "\\end{table}\n");
            }
          fclose (f);

        }
      else
        {
          printf ("\nStarting another calculation...\n");
          std::unique_ptr<mesh> main_mesh = fillers::fill_mesh_by_arguments (argc, argv);

          timed_discrete_function tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get (), "H");
          timed_discrete_function tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V1");
          timed_discrete_function tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V2");
          trio essential (tdfH, tdfV1, tdfV2);
          fillers::fill_initial_info (essential);
          print_parameters (essential);

          if (KNOWN_FUNC)
            {
              timed_discrete_function real_tdfH (main_mesh->m_H_grid.get (), main_mesh->m_scale.get (), "H");
              timed_discrete_function real_tdfV1 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V1");
              timed_discrete_function real_tdfV2 (main_mesh->m_V_grid.get (), main_mesh->m_scale.get (), "V2");
              trio real (real_tdfH, real_tdfV1, real_tdfV2);
              fillers::fill_real_info (real);

              time_loop (essential, &real);

              print_to_gnuplot (tdfH, "H");
              print_to_gnuplot (tdfV1, "V1");
              print_to_gnuplot (tdfV2, "V2");
              print_velocity_to_gnuplot (tdfV1, tdfV2, "V1V2");

              print_to_gnuplot (real_tdfH, "realH");
              print_to_gnuplot (real_tdfV1, "realV1");
              print_to_gnuplot (real_tdfV2, "realV2");
              print_velocity_to_gnuplot (real_tdfV1, real_tdfV2, "realV1V2");

              print_residuals (tdfH, real_tdfH, "H");
              print_residuals (tdfV1, real_tdfV1, "V1");
              print_residuals (tdfV2, real_tdfV2, "V2");
            }
          else
            {
              time_loop (essential, nullptr);

              print_to_gnuplot (tdfH, "H");
              print_to_gnuplot (tdfV1, "V1");
              print_to_gnuplot (tdfV2, "V2");
              print_velocity_to_gnuplot (tdfV1, tdfV2, "V1V2");
            }

          return 0;

        }
    }
  return 0;

}
