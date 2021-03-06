#include "misc.h"

void assert (bool check, std::string message)
{
  if (!check)
    {
      fprintf (stderr, "WARNING: %s\n", message.c_str());
      throw std::runtime_error ("ASSERT: " + message);
    }
}

int fuzzycmp (double a, double b, double eps)
{
  if (fabs (a - b) < eps)
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
