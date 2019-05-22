#pragma once
#include "misc.h"

struct mesh
{
  std::unique_ptr<grid> m_H_grid;
  std::unique_ptr<grid> m_V_grid;
  std::unique_ptr<scale> m_scale;
};

struct trio
{
  trio (const trio &) = delete;
  trio (timed_discrete_function &tdfH, timed_discrete_function &tdfV1, timed_discrete_function &tdfV2)
    : m_tdfH (tdfH), m_tdfV1 (tdfV1), m_tdfV2 (tdfV2) {}
  timed_discrete_function &m_tdfH;
  timed_discrete_function &m_tdfV1;
  timed_discrete_function &m_tdfV2;
};

namespace fillers
{
void fill_initial_info (trio &essential);
void fill_real_info (trio &real);
std::unique_ptr<mesh> fill_mesh_by_arguments (int argc, char **argv);
}
