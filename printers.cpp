#include "printers.h"
#include "discrete_function.h"
#include "fillers_misc.h"
#include "grid.h"

#define SCALE 1

void print_to_gnuplot (discrete_function &df, std::string name, double t)
{
  FILE *out = fopen (name.c_str (), "w");

  if (!out)
    {
      printf ("Cannot create output files\n");
      return;
    }

  df.do_for_each ([out, t] (index ij, point xy, discrete_function &self) {
    double value = self.get_value (ij);
    fprintf (out, "%f %f %f #%f\n", xy.first, xy.second, value, t);
  });

  fclose (out);
}

void print_velocity_to_gnuplot (discrete_function &dfV1, discrete_function &dfV2, std::string name, double t)
{
  FILE *out = fopen (name.c_str (), "w");

  if (!out)
    {
      printf ("Cannot create output files\n");
      return;
    }

  dfV1.do_for_each ([&] (index ij, point xy, discrete_function &) {
    double valueV1 = dfV1.get_value (ij) * SCALE;
    double valueV2 = dfV2.get_value (ij) * SCALE;
    fprintf (out, "%f %f %f %f #%f\n", xy.first, xy.second, valueV1, valueV2, t);
  });

  fclose (out);
}

void print_to_gnuplot (timed_discrete_function &tdf, std::string name)
{
  tdf.do_for_each ([name] (int k, double t, timed_discrete_function &self) {
      std::string filename = name + "/" + std::to_string (k);
      discrete_function &df = self.get_cut (k);
      print_to_gnuplot (df, filename, t);
    });
}

void print_velocity_to_gnuplot (timed_discrete_function &tdfV1, timed_discrete_function &tdfV2, std::string name)
{
  tdfV1.do_for_each ([&] (int k, double t, timed_discrete_function &) {
      std::string filename = name + "/" + std::to_string (k);
      discrete_function &dfV1 = tdfV1.get_cut (k);
      discrete_function &dfV2 = tdfV2.get_cut (k);
      print_velocity_to_gnuplot (dfV1, dfV2, filename, t);
    });
}

void print_residuals (timed_discrete_function &tdf, timed_discrete_function &tdf_real, std::string name)
{
  double val = tdf.residual_C (tdf_real);
  printf ("Residual for %s is: %e\n", name.c_str (), val);
}

void print_parameters (const trio &essential)
{
  const timed_discrete_function &V1 = essential.m_tdfV1;
  const scale *sc = V1.get_scale ();
  scale_parameters sc_p = sc->get_parameters ();
  int t_point_count = sc_p.m_t_point_count;
  const grid *gr = V1.get_cut (0).get_grid ();
  grid_parameters gr_p = gr->get_parameters ();
  int x_point_count = gr_p.m_x_point_count;
  int y_point_count = gr_p.m_y_point_count;
  double miu = MIU;
  printf ("MIU=%.3f,N=%d,M1=%d,M2=%d\n", miu, t_point_count, x_point_count, y_point_count);
}
