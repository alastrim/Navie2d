#include "matvec.h"

int solve_system (std::map<unsigned int, double> &A, std::vector<double> &B, std::vector<double> &X, std::vector<double> &real)
{
  static int print = 1;
  matrix LA;
  vector LB;
  vector LX;
  vector Lreal;

  LA.set (A, toi (B.size ()));
  LB.set (B);
  LX.set (X);
  Lreal.set (real);


  SetRTCAccuracy (EPS_FOR_SOLVING);
  BiCGIter (LA.get_as_laspack (), LX.get_as_laspack (), LB.get_as_laspack (), MAXITER, JacobiPrecond, 0);
//  CGNIter (LA.get_as_laspack (), LX.get_as_laspack (), LB.get_as_laspack (), MAXITER, nullptr, 0);

  if (DEBUG && print > 0)
    {
//      printf ("Matrix:\n");
//      LA.print ();
//      printf ("Vector:\n");
//      LB.print ();
      printf ("System:\n");
      print_system (LA, LB);
      printf ("Result:\n");
      LX.print ();
      printf ("Real:\n");
      Lreal.print ();
      printf ("\n");
      print--;
    }

  X = LX.get ();
  return 0;
}

// matrix----------------------------------------------------------------------

matrix::matrix () : m_name ("Empty matrix"), m_size (0), m_non_zero_count (0), m_laspack_pointer (nullptr)
{
}

matrix::~matrix ()
{
  drop_laspack_pointer ();
}

void matrix::set (std::map<unsigned int, double> &src, int full_size)
{
  m_name = std::string ("Non-empty matrix");
  m_container = src;
  m_size = full_size;
  m_non_zero_count = m_size + 1;
  m_non_zero_count = MAX_NON_ZERO;
}

QMatrix *matrix::get_as_laspack ()
{
  update_to_laspack ();
  return m_laspack_pointer;
}

int matrix::size ()
{
  return m_size;
}

void matrix::print ()
{
  for (int i = 0; i < m_size; i++)
    {
      for (int j = 0; j < m_size; j++)
        {
          double value = 0;
          if (m_container.find (tou (i * m_size + j)) != m_container.end ())
            value = m_container[tou (i * m_size + j)];
          printf ("%9.6f ", value);
        }
      printf ("\n");
    }
}

void matrix::drop_laspack_pointer ()
{
  if (m_laspack_pointer)
    Q_Destr (m_laspack_pointer);

  m_laspack_pointer = nullptr;
}

void matrix::update_to_laspack ()
{
  drop_laspack_pointer ();
  m_laspack_pointer = new QMatrix ();
  Q_Constr (m_laspack_pointer, const_cast<char *> (m_name.c_str ()), tou (m_size), False, Rowws, Normal, True);

  for (int li = 1; li <= m_size; li++)
    Q_SetLen (m_laspack_pointer, tou (li), tou (m_non_zero_count));

  std::vector<int> non_zero_counts (tou (m_size), 0);
  for (const std::pair<unsigned int, double> ind_and_val : m_container)
    {
      int i = toi (ind_and_val.first) / m_size;
      int j = toi (ind_and_val.first) - i * m_size;

      int li = i + 1;
      int lj = j + 1;

      double value = ind_and_val.second;
      Q_SetEntry (m_laspack_pointer, tou (li), tou (non_zero_counts[tou (i)]++), tou (lj), value);
    }
}

// vector----------------------------------------------------------------------

vector::vector () : m_name ("Empty vector"), m_size (0), m_laspack_pointer (nullptr)
{
}

vector::~vector ()
{
  drop_laspack_pointer ();
}

void vector::set (std::vector<double> &src)
{
  m_name = "Vector from std";
  m_container = src;
  m_size = toi (m_container.size ());
}

std::vector<double> &vector::get ()
{
  update_from_laspack ();
  return m_container;
}

Vector *vector::get_as_laspack ()
{
  update_to_laspack ();
  return m_laspack_pointer;
}

int vector::size ()
{
  return m_size;
}

void vector::print ()
{
  update_from_laspack ();

  for (int i = 0; i < m_size; i++)
      printf ("%9.6f ", m_container[tou (i)]);
  printf ("\n");
}

void vector::drop_laspack_pointer ()
{
  if (m_laspack_pointer)
    V_Destr (m_laspack_pointer);
  m_laspack_pointer = nullptr;
}

void vector::update_to_laspack ()
{
  drop_laspack_pointer ();
  m_laspack_pointer = new Vector ();
  V_Constr (m_laspack_pointer, const_cast<char *> (m_name.c_str ()), tou (m_size), Normal, True);

  for (int i = 0; i < m_size; i++)
    {
      int li = i + 1;
      double value = m_container[tou (i)];
      V_SetCmp (m_laspack_pointer, tou (li), value);
    }
}

void vector::update_from_laspack ()
{
  if (!m_laspack_pointer)
    return;
  for (int i = 0; i < m_size; i++)
    {
      int li = i + 1;
      m_container[tou (i)] = V_GetCmp (m_laspack_pointer, tou (li));
    }
}

// ----------------------------------------------------------------------


void print_system (matrix &a, vector &b)
{
  b.update_from_laspack ();
  for (int i = 0; i < a.m_size; i++)
    {
      for (int j = 0; j < a.m_size; j++)
        printf ("%9.6f ", a.m_container[tou (i * a.m_size + j)]);
      printf (" %9.6f\n", b.m_container[tou (i)]);
    }
}

