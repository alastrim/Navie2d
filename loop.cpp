#include "loop.h"
#include "discrete_function.h"
#include "grid.h"
#include "matvec.h"
#include "fillers.h"

void time_loop (trio &essential)
{
  essential.m_tdfH.do_for_each ([&essential] (int k, double, timed_discrete_function &)
  {
      {
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> &X = essential.m_tdfH.get_cut (k + 1).get_raw_vector ();

        fillers::fill_first (k, A, B, essential);

        solve_system (A, B, X);
      }
      // Solving second system
      {
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> &X = essential.m_tdfV1.get_cut (k + 1).get_raw_vector ();

        fillers::fill_second (k, A, B, essential);

        solve_system (A, B, X);
      }
      // Solving third system
      {
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> &X = essential.m_tdfV2.get_cut (k + 1).get_raw_vector ();

        fillers::fill_third (k, A, B, essential);

        solve_system (A, B, X);
      }
    });
}
