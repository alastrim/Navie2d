#include "loop.h"
#include "discrete_function.h"
#include "grid.h"
#include "matvec.h"
#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"

void time_loop (trio &essential, trio &real)
{
  essential.m_tdfH.do_for_each ([&essential, &real] (int k, double, timed_discrete_function &)
  {
      if (k == essential.m_tdfH.get_scale ()->get_parameters ().m_t_step_count)
        return;

      {
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> &X = essential.m_tdfH.get_cut (k + 1).get_raw_vector ();
        std::vector<double> &realvect = real.m_tdfH.get_cut (k + 1).get_raw_vector ();

        fillers::fill_first (k, A, B, essential);

        solve_system (A, B, X, realvect);
      }
      // Solving second system
      {
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> &X = essential.m_tdfV1.get_cut (k + 1).get_raw_vector ();
        std::vector<double> &realvect = real.m_tdfV1.get_cut (k + 1).get_raw_vector ();

        fillers::fill_second (k, A, B, essential);

        solve_system (A, B, X, realvect);
      }
      // Solving third system
      {
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> &X = essential.m_tdfV2.get_cut (k + 1).get_raw_vector ();
        std::vector<double> &realvect = real.m_tdfV2.get_cut (k + 1).get_raw_vector ();

        fillers::fill_third (k, A, B, essential);

        solve_system (A, B, X, realvect);
      }
    });
}
