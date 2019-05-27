#pragma once
#include "misc.h"

void print_to_gnuplot (discrete_function &df, std::string name, double t);
void print_to_gnuplot (timed_discrete_function &tdf, std::string name);
void print_residuals (std::string filename, timed_discrete_function &tdf, timed_discrete_function &tdf_real, std::string name);
std::string get_parameters_string (const trio &essential);
