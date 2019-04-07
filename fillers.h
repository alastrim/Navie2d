#ifndef FILLERS_H
#define FILLERS_H

#include "misc.h"

//void fill_first (std::vector<double> &A, std::vector<double> &B, discrete_function &H, discrete_function &V, int n, double h, double tau, discrete_function &f, discrete_function &f_0);
//void fill_second (std::vector<double> &A, std::vector<double> &B, discrete_function &H, discrete_function &V, int n, double h, double tau, discrete_function &f, discrete_function &f_0);

struct mesh
{
  std::unique_ptr<grid> m_H_grid;
  std::unique_ptr<grid> m_V_grid;
  std::unique_ptr<scale> m_scale;
};

struct trio
{
  trio (timed_discrete_function &tdfH, timed_discrete_function &tdfV1, timed_discrete_function &tdfV2)
    : m_tdfH (tdfH), m_tdfV1 (tdfV1), m_tdfV2 (tdfV2) {}
  timed_discrete_function &m_tdfH;
  timed_discrete_function &m_tdfV1;
  timed_discrete_function &m_tdfV2;
};

namespace fillers
{
void fill_first (int k, std::vector<double> &A, std::vector<double> &B, trio &essential);
void fill_second (int k, std::vector<double> &A, std::vector<double> &B, trio &essential);
void fill_third (int k, std::vector<double> &A, std::vector<double> &B, trio &essential);

void fill_initial_info (trio essential);
std::unique_ptr<mesh> fill_mesh_by_arguments (int argc, char **argv);
}

#endif // FILLERS_H
