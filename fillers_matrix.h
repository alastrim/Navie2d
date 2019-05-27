#pragma once
#include "misc.h"

namespace fillers
{
void fill_first (int k, std::unordered_map<unsigned int, double> &A, std::vector<double> &B, trio &essential);
void fill_second (int k, std::unordered_map<unsigned int, double> &A, std::vector<double> &B, trio &essential);
void fill_third (int k, std::unordered_map<unsigned int, double> &A, std::vector<double> &B, trio &essential);
}
