#ifndef MISC_H
#define MISC_H

extern "C"
{
#include <laspack/version.h>
#include <laspack/mlsolv.h>
#include <laspack/elcmp.h>
#include <laspack/copyrght.h>
#include <laspack/rtc.h>
#include <laspack/eigenval.h>
#include <laspack/precond.h>
#include <laspack/matrix.h>
#include <laspack/factor.h>
#include <laspack/lastypes.h>
#include <laspack/errhandl.h>
#include <laspack/qmatrix.h>
#include <laspack/itersolv.h>
#include <laspack/vector.h>
#include <laspack/operats.h>
#include <xc/getopts.h>
#include <xc/xtypes.h>
#include <xc/version.h>
#include <xc/xstring.h>
}
#include <string>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <functional>
#include <memory>

#define MIN_FOR_DIVISION 1e-16
#define MIN_FOR_COMPARISON 1e-16
#define EPS_FOR_SOLVING 1e-5
#define MAXITER 100
#define MAX_NON_ZERO 3
#define ELEMS_ON_SCREEN 100
#define GAMMA 1.4
#define MIU 0.001
#define DEBUG 1

void assert (bool check, std::string message);
int toi (size_t src);
unsigned int tou (int src);
int fuzzycmp (double a, double b = 0.0);
unsigned int sbsc (int step_count);

class grid;
class scale;
class discrete_function;
class timed_discrete_function;
struct trio;

typedef std::pair<double, double> point;
typedef std::pair<int, int> index;
typedef std::function<double (point xy)> continuous_function;
typedef std::function<double (double t, point xy)> timed_continuous_function;
typedef std::function<void (index ij, point xy, discrete_function &self)> discrete_foreach_function;
typedef std::function<void (int k, double t, timed_discrete_function &self)> timed_discrete_foreach_function;

#endif // MISC_H
