#include "fillers.h"
#include "misc.h"
#include "matvec.h"
#include "discrete_function.h"
#include "grid.h"

namespace fillers
{
void fill_first (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
    std::vector<double> &H_raw = essential.m_tdfH.get_cut (k).get_raw_vector ();
    std::vector<double> &V1_raw = essential.m_tdfV1.get_cut (k).get_raw_vector ();
    std::vector<double> &V2_raw = essential.m_tdfV2.get_cut (k).get_raw_vector ();
    unsigned int S = static_cast<unsigned int> (H_raw.size ());

    A = std::vector<double> (S * S, 0);
    B = std::vector<double> (S, 0);

    for (unsigned int s = 0; s < S; s++)
      {
        A[s * S + s] = 1.5;
        B[s] = H_raw[s] + H_raw[s] * 0 + V1_raw[s] * 0 + V2_raw[s] * 0;
      }
}

void fill_second (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
    std::vector<double> &H_raw = essential.m_tdfH.get_cut (k).get_raw_vector ();
    std::vector<double> &V1_raw = essential.m_tdfV1.get_cut (k).get_raw_vector ();
    std::vector<double> &V2_raw = essential.m_tdfV2.get_cut (k).get_raw_vector ();
    unsigned int S = static_cast<unsigned int> (V1_raw.size ());

    A = std::vector<double> (S * S, 0);
    B = std::vector<double> (S, 0);

    for (unsigned int s = 0; s < S; s++)
      {
        A[s * S + s] = 1.5;
        B[s] = V1_raw[s] + H_raw[s] * 0 + V1_raw[s] * 0 + V2_raw[s] * 0;
      }
}

void fill_third (int k, std::vector<double> &A, std::vector<double> &B, trio &essential)
{
    std::vector<double> &H_raw = essential.m_tdfH.get_cut (k).get_raw_vector ();
    std::vector<double> &V1_raw = essential.m_tdfV1.get_cut (k).get_raw_vector ();
    std::vector<double> &V2_raw = essential.m_tdfV2.get_cut (k).get_raw_vector ();
    unsigned int S = static_cast<unsigned int> (V2_raw.size ());

    A = std::vector<double> (S * S, 0);
    B = std::vector<double> (S, 0);

    for (unsigned int s = 0; s < S; s++)
      {
        A[s * S + s] = 1.5;
        B[s] = V2_raw[s] + H_raw[s] * 0 + V1_raw[s] * 0 + V2_raw[s] * 0;
      }
}

void fill_initial_info (trio essential)
{
  discrete_function &H_initial_cut = essential.m_tdfH.get_cut (0);
  discrete_function &V1_initial_cut = essential.m_tdfV1.get_cut (0);
  discrete_function &V2_initial_cut = essential.m_tdfV2.get_cut (0);

  continuous_function H_initial_filler = [] (point xy) { double x = xy.first, y = xy.second; return x * y; };
  continuous_function V1_initial_filler = [] (point xy) { double x = xy.first, y = xy.second; return x + y; };
  continuous_function V2_initial_filler = [] (point xy) { double x = xy.first, y = xy.second; return x - y; };

  H_initial_cut.fill (H_initial_filler);
  V1_initial_cut.fill (V1_initial_filler);
  V2_initial_cut.fill (V2_initial_filler);

  discrete_foreach_function zero_setter = [&] (index ij, point, discrete_function &self) { self.set_value (ij, 0); };
  timed_discrete_foreach_function timed_zero_setter = [&] (int k, double, timed_discrete_function &self) { self.get_cut (k).do_for_edge (zero_setter); };

  essential.m_tdfH.do_for_each (timed_zero_setter);
  essential.m_tdfV1.do_for_each (timed_zero_setter);
  essential.m_tdfV2.do_for_each (timed_zero_setter);
}

std::unique_ptr<mesh> fill_mesh_by_arguments (int argc, char **argv)
{
  int t_step_count, x_step_count, y_step_count, iT, iX, iY;
  if (argc < 7 || !(iT = atoi (argv[1])) || !(t_step_count = atoi(argv[2])) || !(iX = atoi (argv[3])) || !(x_step_count = atoi(argv[4])) || !(iY = atoi (argv[5])) || !(y_step_count = atoi(argv[6])))
    {
      printf ("Usage: ./main.exe <T> <t_step_count> <X> <x_step_count> <Y> <y_step_count>\n");
      printf ("Bad arguments given, using default...");
      t_step_count = 98;
      x_step_count = 55;
      y_step_count = 76;
      iT = 9;
      iX = 2;
      iY = 3;
    }
  double T = iT, X = iX, Y = iY;

  std::unique_ptr<mesh> result = std::make_unique<mesh> ();

  result->m_H_grid = std::make_unique<grid> (0, 0, X, Y, x_step_count, y_step_count);
  result->m_V_grid = std::make_unique<grid> (0, 0, X, Y, x_step_count, y_step_count);
  result->m_scale = std::make_unique<scale> (0, T, t_step_count);

  return result;
}
}

//void fill_first (std::vector<double> &A, std::vector<double> &B,
//                 discrete_function &H, discrete_function &V,
//                 int n, double h, double tau, discrete_function &/*f*/,
//                 discrete_function &f_0)
//{
//  std::vector<double> &H_cut = H.cut (n);
//  std::vector<double> &V_cut = V.cut (n);
//  unsigned int M = static_cast<unsigned int> (H_cut.size ());

//  A = std::vector<double> (M * M, 0);
//  B = std::vector<double> (M, 0);

//  for (unsigned int m = 0; m < M; m++)
//    {
//      if (m > 0)
//        A[m * M + m - 1] = -(V_cut[m] + fabs (V_cut[m])) / (2.0 * h);

//      A[m * M + m] = 1.0 / tau + (V_cut[m + 1] + fabs (V_cut[m + 1])
//          - V_cut[m] + fabs (V_cut[m])) / (2.0 * h);

//      if (m < M - 1)
//        A[m * M + m + 1] = (V_cut[m + 1] - fabs (V_cut[m + 1])) / (2.0 * h);

//      B[m] = H_cut[m] * (1.0 / tau) + f_0.val (n, toi (m));
//    }
//}

//void fill_second (std::vector<double> &A, std::vector<double> &B,
//                  discrete_function &H, discrete_function &V,
//                  int n, double h, double tau, discrete_function &f,
//                  discrete_function &/*f_0*/)
//{
//  std::vector<double> &H_cut = H.cut (n + 1);
//  std::vector<double> &V_cut = V.cut (n);
//  unsigned int M = static_cast<unsigned int> (H_cut.size ());

//  A = std::vector<double> ((M + 1) * (M + 1), 0);
//  B = std::vector<double> (M + 1, 0);

//  for (unsigned int m = 0; m < M + 1; m++)
//    {
//      double H_s_;
//      if (m == 0 || m == M || !fuzzycmp ((H_s_ = (H_cut[m] + H_cut[m - 1]) / 2)))
//        {
//          A[m * (M + 1) + m] = 1;
//          B[m] = 0;
//          continue;
//        }

//      A[m * (M + 1) + m - 1] = H_s_ * -(V_cut[m] + fabs (V_cut[m]))
//                               / (2 * h) - MIU / (h * h);

//      A[m * (M + 1) + m] = H_s_ * (1.0 / tau + (V_cut[m] + fabs (V_cut[m]))
//                                   / (2 * h) - (V_cut[m] - fabs (V_cut[m]))
//                                   / (2 * h)) + (2 * MIU) / (h * h);

//      A[m * (M + 1) + m + 1] = H_s_ * (V_cut[m] - fabs (V_cut[m]))
//                               / (2 * h) - MIU / (h * h);

//      B[m] = H_s_ * (V_cut[m] / tau - (GAMMA / (GAMMA - 1)) *
//                     ((pow (H_cut[m], GAMMA - 1) - pow (H_cut[m - 1], GAMMA - 1))
//                     / h) + f.val (n + 1, toi (m)));
//    }
//}
