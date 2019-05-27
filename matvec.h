#pragma once
#include "misc.h"
class matrix;
class vector;

// Interfaces from C++ containers to laspack matrix and vector

int solve_system (std::unordered_map<unsigned int, double> &A, std::vector<double> &B, std::vector<double> &X, std::vector<double> &real);

class matrix
{
friend void print_system (matrix &a, vector &b);
public:
  matrix ();
  ~matrix ();

  void set (std::unordered_map<unsigned int, double> &src, int full_size);
  // std::vector<double> &get ();

  QMatrix *get_as_laspack ();
  int size ();
  void print ();

private:
  void drop_laspack_pointer ();
  void update_to_laspack ();
  // void update_from_laspack ();

  std::string m_name;
  int m_size;
  int m_non_zero_count;
  std::unordered_map<unsigned int, double> m_container;
  QMatrix *m_laspack_pointer;
};

class vector
{
friend void print_system (matrix &a, vector &b);
public:
  vector ();
  ~vector ();

  void set (std::vector<double> &src);
  std::vector<double> &get ();

  Vector *get_as_laspack ();
  int size ();
  void print ();

private:
  void drop_laspack_pointer ();
  void update_to_laspack ();
  void update_from_laspack ();

  std::string m_name;
  int m_size = 0;
  std::vector<double> m_container;
  Vector *m_laspack_pointer;
};

void print_system (matrix &a, vector &b);
