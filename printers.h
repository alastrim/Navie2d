#ifndef PRINTERS_H
#define PRINTERS_H

#include "misc.h"

void print_to_gnuplot (discrete_function &df, std::string name);
void print_to_gnuplot (timed_discrete_function &tdf, std::string name);

#endif // PRINTERS_H
