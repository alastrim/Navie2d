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

  df.do_for_each ([out] (index ij, point xy, discrete_function &self) {
    double value = self.get_value (ij);
    fprintf (out, "%f %f %f\n", xy.first, xy.second, value);
  });

  fclose (out);
}

void print_to_gnuplot (timed_discrete_function &tdf, std::string name)
{
  tdf.do_for_each ([name] (int k, double, timed_discrete_function &self) {
      std::string filename = name + "/" + std::to_string (k);
      discrete_function &df = self.get_cut (k);
      print_to_gnuplot (df, filename);
    });
}
