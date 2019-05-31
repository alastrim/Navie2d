#include "discrete_function.h"
#include "grid.h"

discrete_function::discrete_function (const grid *grid, std::string name) : m_grid (grid)
{
  m_i_size = tou (m_grid->get_parameters ().m_x_point_count);
  m_j_size = tou (m_grid->get_parameters ().m_y_point_count);

  m_data.resize (m_i_size * m_j_size);
  fill ([] (point) { return 0; });
  m_name = name;
}

void discrete_function::fill (continuous_function cf)
{
//  printf ("%s\n", m_name.c_str ());
//  for (int i = 0; i < toi (m_i_size); i++)
//    {
//      for (int j = 0; j < toi (m_j_size); j++)
//        {
//          point_type type = m_grid->get_type ({i, j});
//          switch (type)
//            {
//          case point_type::inner:
//          {
//            printf ("i");
//            break;
//          }
//          case point_type::edge:
//          {
//            printf ("e");
//            break;
//          }
//          case point_type::outer:
//          {
//            printf ("o");
//            break;
//          }
//          case point_type::INVALID:
//          {
//            assert (false, "Sanity");
//          }
//            }
//        }
//      printf ("\n");
//    }
//  printf ("\n");

  do_for_each ([&cf] (index ij, point xy, discrete_function &self){
    double value = cf (xy);
    if (self.m_grid->get_full_type (ij) != point_type::outer)
      self.set_value (ij, value);
  });
}

double discrete_function::residual (const discrete_function &real)
{
  double max = 0;
  do_for_each ([&max, &real] (index ij, point, discrete_function &self){
      double diff = fabs (self.get_value (ij) - real.get_value (ij));
      if (diff > max)
        max = diff;
  });
  return max;
}

void discrete_function::set_value (index ij, double value)
{
  assert (m_grid->get_full_type (ij) != point_type::outer, "Bad argument for function set value");

  m_data[tou (ij.first) * m_j_size + tou (ij.second)] = value;
}

double discrete_function::get_value (index ij) const
{
  if (m_grid->get_full_type (ij) == point_type::outer)
    return 0;

  return m_data[tou (ij.first) * m_j_size + tou (ij.second)];
}

double discrete_function::left (index ij) const
{
  assert (m_name == "H", "Do not call left/right for anything other than H");
  return (get_value (ij) + get_value ({ij.first, ij.second - 1})) / 2.;
}

double discrete_function::right (index ij) const
{
  assert (m_name == "H", "Do not call left/right for anything other than H");
  return (get_value (ij) + get_value ({ij.first - 1, ij.second})) / 2.;
}

double discrete_function::tilda (index ij) const
{
  if (m_name == "V1")
    {
      return (get_value (ij) + get_value ({ij.first, ij.second + 1})) / 2.;
    }
  if (m_name == "V2")
    {
      return (get_value (ij) + get_value ({ij.first + 1, ij.second})) / 2.;
    }
  assert (false, "Do not call tilda for anything other than V1 and V2");
  return 0;
}

void discrete_function::do_for_each (discrete_foreach_function dff)
{
  for (int i = 0; i < toi (m_i_size); i++)
    {
      for (int j = 0; j < toi (m_j_size); j++)
        {
          index ij = {i, j};
//          if (m_grid->get_type (ij) == point_type::outer)
//            continue;
          point xy = m_grid->get_point (ij);
          dff (ij, xy, *this);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

timed_discrete_function::timed_discrete_function (const grid *grid, const scale *scale, std::string name)
  : m_name (name), m_grid (grid), m_scale (scale)
{
  m_k_size = tou (m_scale->get_parameters ().m_t_point_count);

  m_data.resize (m_k_size);
  fill ([] (double, point) { return 0; });
}

void timed_discrete_function::fill (timed_continuous_function tcf)
{
  do_for_each ([&tcf] (int k, double t, timed_discrete_function &self) {
      continuous_function cf = [&] (point xy) { return tcf (t, xy); };
      std::unique_ptr<discrete_function> df = std::make_unique<discrete_function> (self.m_grid, self.m_name);
      df->fill (cf);
      self.set_cut (k, std::move (df));
    });
}

double timed_discrete_function::residual (const timed_discrete_function &real)
{
  double max = 0;
  do_for_each ([&max, &real] (int k, double, timed_discrete_function &self) {
      double diff = self.get_cut (k).residual (real.get_cut (k));
      if (diff > max)
        max = diff;
    });
  return max;
}

void timed_discrete_function::set_cut (int k, std::unique_ptr<discrete_function> df)
{
  assert (k >= 0 && tou (k) < m_k_size, "Bad argument for timed function cut");

  std::unique_ptr<discrete_function> &cut = m_data[tou (k)];
  cut = std::move (df);
}

discrete_function &timed_discrete_function::get_cut (int k)
{
  assert (k >= 0 && tou (k) < m_k_size, "Bad argument for timed function cut");

  return *m_data[tou (k)];
}

const discrete_function &timed_discrete_function::get_cut (int k) const
{
  assert (k >= 0 && tou (k) < m_k_size, "Bad argument for timed function cut");

  return *m_data[tou (k)];
}

void timed_discrete_function::do_for_each (timed_discrete_foreach_function tdff)
{
  for (int k = 0; k < toi (m_k_size); k++)
    {
      double t = m_scale->get_time (k);
      tdff (k, t, *this);
    }
}

