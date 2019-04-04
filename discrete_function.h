#ifndef DICRETE_FUNCTION_H
#define DICRETE_FUNCTION_H

#include "misc.h"

class discrete_function
{
public:
  discrete_function (const grid *grid);
  void fill (continuous_function cf);
  void set_value (index ij, double value);
  double get_value (index ij);
  const grid *get_grid () { return m_grid; }
  void do_for_each (discrete_foreach_function dff);
private:
  const grid *m_grid;
  std::vector<double> m_data;
  unsigned int m_i_size;
  unsigned int m_j_size;
};

class timed_discrete_function
{
public:
  timed_discrete_function (const grid *grid, const scale *scale);
  void fill (timed_continuous_function tcf);
  void set_cut (int k, std::unique_ptr<discrete_function> df);
  discrete_function &get_cut (int k);
  const grid *get_grid () { return m_grid; }
  const scale *get_scale () { return m_scale; }
  void do_for_each (timed_discrete_foreach_function tdff);
private:
  const grid *m_grid;
  const scale *m_scale;
  std::vector<std::unique_ptr<discrete_function>> m_data;
  unsigned int m_k_size;
};

#endif // DICRETE_FUNCTION_H
