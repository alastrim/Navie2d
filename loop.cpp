#include "loop.h"
#include "discrete_function.h"
#include "grid.h"
#include "matvec.h"
#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"

static void set_V1_edge_condition (discrete_function &df)
{
  if (KNOWN_FUNC)
    return;

  df.do_for_each ([&] (index ij, point, discrete_function &)
  {
    if (df.get_grid ()->get_full_type (ij) == point_type::outer)
      return;

    if (ij.first == 0)
      df.set_value (ij, OMEGA);
  });
}

static void set_H_edge_condition (discrete_function &df)
{
  if (KNOWN_FUNC)
    return;

  df.do_for_each ([&] (index ij, point, discrete_function &)
  {
    if (df.get_grid ()->get_full_type (ij) == point_type::outer)
      return;

    if (ij.first == 0)
      df.set_value (ij, RHO_GAMMA);
  });
}

void time_loop (trio &essential, trio *real)
{
  essential.m_tdfH.do_for_each ([&essential, real] (int k, double, timed_discrete_function &)
  {
      printf ("Looping step %d of %d\n", k, essential.m_tdfH.get_scale ()->get_parameters ().m_t_point_count - 1);
      if (k == essential.m_tdfH.get_scale ()->get_parameters ().m_t_point_count - 1)
        return;

      set_V1_edge_condition (essential.m_tdfV1.get_cut (k));
      set_H_edge_condition (essential.m_tdfH.get_cut (k));

      {
        std::vector<double> &X = essential.m_tdfH.get_cut (k + 1).get_raw_vector ();
        std::map<unsigned int, double> A;
        std::vector<double> B = std::vector<double> (X.size (), 0);

        fillers::fill_first (k, A, B, essential);

        if (real)
          {
            std::vector<double> &realvect = real->m_tdfH.get_cut (k + 1).get_raw_vector ();
            solve_system (A, B, X, &realvect);
          }
        else
          solve_system (A, B, X, nullptr);
      }

      // Solving second system
      {
        std::vector<double> &X = essential.m_tdfV1.get_cut (k + 1).get_raw_vector ();
        std::map<unsigned int, double> A;
        std::vector<double> B = std::vector<double> (X.size (), 0);

        fillers::fill_second (k, A, B, essential);

        if (real)
          {
            std::vector<double> &realvect = real->m_tdfV2.get_cut (k + 1).get_raw_vector ();
            solve_system (A, B, X, &realvect);
          }
        else
          solve_system (A, B, X, nullptr);
      }

      // Solving third system
      {
        std::vector<double> &X = essential.m_tdfV2.get_cut (k + 1).get_raw_vector ();
        std::map<unsigned int, double> A;
        std::vector<double> B = std::vector<double> (X.size (), 0);

        fillers::fill_third (k, A, B, essential);

        if (real)
          {
            std::vector<double> &realvect = real->m_tdfV2.get_cut (k + 1).get_raw_vector ();
            solve_system (A, B, X, &realvect);
          }
        else
          solve_system (A, B, X, nullptr);
      }

      set_V1_edge_condition (essential.m_tdfV1.get_cut (k + 1));
      set_H_edge_condition (essential.m_tdfH.get_cut (k + 1));
    });
}
