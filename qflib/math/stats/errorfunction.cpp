/**
@file  errorfunction.cpp
@brief Implementation of the inverse complementary error function
*/

#include <qflib/math/stats/errorfunction.hpp>
#include <cmath>

BEGIN_NAMESPACE(qf)

double ErrorFunction::inverfc(double p)
{
  // Inverse of complementary error function.
  // Returns x such that erfc(x) = p for argument p between 0 and 2.
  if (p >= 2.0)
    return -100.;
  if (p <= 0.0)
    return 100.;
  double pp = (p < 1.0) ? p : 2. - p;
  double t = std::sqrt(-2. * std::log(pp / 2.));
  // Initial guess from rational approximation:
  double x = -0.70711 * ((2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t);
  // Halley's method refinement: cubic convergence per iteration.
  // erfc'(x) = -2/sqrt(pi) * exp(-x^2), so 1.12837916709551257 = 2/sqrt(pi).
  for (int j = 0; j < 4; j++) {
    double err = std::erfc(x) - pp;
    x += err / (1.12837916709551257 * std::exp(-x * x) - x * err);
    x = x < 0 ? 0 : x;
  }
  return (p < 1.0 ? x : -x);
}

END_NAMESPACE(qf)
