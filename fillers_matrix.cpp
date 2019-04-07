#include "fillers_misc.h"
#include "fillers_matrix.h"
#include "fillers_functions.h"
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
        A[s * S + s] = 1.1;
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
        A[s * S + s] = 1.1;
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
        A[s * S + s] = 1.1;
        B[s] = V2_raw[s] + H_raw[s] * 0 + V1_raw[s] * 0 + V2_raw[s] * 0;
      }
}
}
