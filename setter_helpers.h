#pragma once
#include "misc.h"

double xmabs (double x);
double xpabs (double x);

class VectorSetter;
class MatrixSetter;

bool process_V_edge (int m1, int m2, double check, MatrixSetter &A_at, VectorSetter &B_at);
bool process_H_edge (int m1, int m2, MatrixSetter &A_at, VectorSetter &B_at);
bool process_V_condition (int m1, int m2, MatrixSetter &A_at, VectorSetter &B_at);

class MatrixSetter
{
friend bool process_V_condition (int m1, int m2, MatrixSetter &A_at, VectorSetter &B_at);
friend bool process_V_edge (int m1, int m2, double check, MatrixSetter &A_at, VectorSetter &B_at);
friend bool process_H_edge (int m1, int m2, MatrixSetter &A_at, VectorSetter &B_at);
public:
  MatrixSetter (std::map<unsigned int, double> &A, const discrete_function &df);
  double & operator () (int m1_base, int m2_base, int m1_mod, int m2_mod);
private:
  std::map<unsigned int, double> &m_A;
  std::vector<int> m_taken;
  const discrete_function &m_df;
  double m_dummy;
};

class VectorSetter
{
public:
  VectorSetter (std::vector<double> &B, const discrete_function &df);
  double & operator () (int m1, int m2);
private:
  std::vector<double> &m_B;
  std::vector<int> m_taken;
  const discrete_function &m_df;
};
