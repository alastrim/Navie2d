﻿#pragma once
#include "misc.h"

class discrete_function
{
public:
  discrete_function (const discrete_function &) = delete;
  discrete_function (const grid *grid, std::string name);
  double residual_C (const discrete_function &real);
  double residual_L2 (const discrete_function &real);
  void fill (continuous_function cf);
  void set_value (index ij, double value);
  double get_value (index ij) const;
  double val (int i, int j) const { return get_value ({i, j}); }
  double tilda (int i, int j) const { return tilda ({i, j}); }
  double left (int i, int j) const { return left ({i, j}); }
  double right (int i, int j) const { return right ({i, j}); }
  double tilda (index ij) const;
  double left (index ij) const;
  double right (index ij) const;
  const grid *get_grid () const { return m_grid; }
  std::string get_name () const { return m_name; }
  void do_for_each (discrete_foreach_function dff);
  std::vector<double> &get_raw_vector () { return m_data; }
private:
  std::string m_name;
  const grid *m_grid;
  std::vector<double> m_data;
  unsigned int m_i_size;
  unsigned int m_j_size;
};

class timed_discrete_function
{
public:
  timed_discrete_function (const timed_discrete_function &) = delete;
  timed_discrete_function (const grid *grid, const scale *scale, std::string name);
  double residual_C (const timed_discrete_function &real);
  double residual_L2 (const timed_discrete_function &real);
  void fill (timed_continuous_function tcf);
  void set_cut (int k, std::unique_ptr<discrete_function> df);
  discrete_function &get_cut (int k);
  const discrete_function &get_cut (int k) const;
  const grid *get_grid () const { return m_grid; }
  const scale *get_scale () const { return m_scale; }
  void do_for_each (timed_discrete_foreach_function tdff);
private:
  std::string m_name;
  const grid *m_grid;
  const scale *m_scale;
  std::vector<std::unique_ptr<discrete_function>> m_data;
  unsigned int m_k_size;
};
