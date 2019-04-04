#include "misc.h"
#include "grid.h"
#include "discrete_function.h"
#include "printers.h"

int main (int argc, char **argv)
{
  int t_step_count, x_step_count, y_step_count, iT, iX, iY;
  if (argc < 7 || !(iT = atoi (argv[1])) || !(t_step_count = atoi(argv[2])) || !(iX = atoi (argv[3])) || !(x_step_count = atoi(argv[4])) || !(iY = atoi (argv[5])) || !(y_step_count = atoi(argv[6])))
    {
      printf ("Usage: ./main.exe <T> <t_step_count> <X> <x_step_count> <Y> <y_step_count>\n");
      return -1;
    }
  double T = iT, X = iX, Y = iY;

  grid H_grid (0, 0, X, Y, x_step_count, y_step_count);
  scale H_scale (0, T, t_step_count);
  timed_discrete_function tdfH (&H_grid, &H_scale);

  tdfH.fill ([] (double t, point xy) {
    double x = xy.first;
    double y = xy.second;
    return sin (t) * cos (x) * exp (y);
  });

  print_to_gnuplot (tdfH, "H");
  return 0;
}
