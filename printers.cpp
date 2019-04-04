#include "printers.h"
#include "discrete_function.h"
#include "grid.h"

void print_to_gnuplot (discrete_function &df, std::string name)
{
  FILE *out = fopen (name.c_str (), "w");

  if (!out)
    {
      printf ("Cannot create output files\n");
      return;
    }

  df.do_for_each ([&] (index ij, point xy) {
    double value = df.get_value (ij);
    fprintf (out, "%f %f %f\n", xy.first, xy.second, value);
  });

  fclose (out);
}

void print_to_gnuplot (timed_discrete_function &tdf, std::string name)
{
  tdf.do_for_each ([&] (int k, double) {
      std::string filename = name + "/" + std::to_string (k);
      discrete_function &df = tdf.get_cut (k);
      print_to_gnuplot (df, filename);
    });
}
