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

namespace fillers
{
void fill_initial_info (timed_discrete_function &tdfH, timed_discrete_function &tdfV1, timed_discrete_function &tdfV2);
std::unique_ptr<mesh> fill_mesh_by_arguments (int argc, char **argv);
}

#endif // FILLERS_H
