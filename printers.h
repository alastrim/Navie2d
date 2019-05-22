#pragma once
#include "misc.h"

void print_to_gnuplot (discrete_function &df, std::string name, double t);
void print_to_gnuplot (timed_discrete_function &tdf, std::string name);
void print_residuals (timed_discrete_function &tdf, timed_discrete_function &tdf_real, std::string name);
