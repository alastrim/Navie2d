#include "misc.h"

void assert (bool check, std::string message)
{
  if (!check)
    fprintf (stderr, "WARNING: %s\n", message.c_str());
}

int fuzzycmp (double a, double b)
{
  if (fabs (a - b) < MIN_FOR_COMPARISON)
    return 0;
  if (a > b)
    return 1;
  return -1;
}

int toi (size_t src)
{
  return static_cast<int> (src);
}

unsigned int tou (int src)
{
  return static_cast<unsigned int> (src);
}

unsigned int sbsc (int step_count)
{
  assert (step_count >= 0, "Bad step count given");
  return tou (step_count + 1);
}

