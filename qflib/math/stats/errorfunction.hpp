/**
@file  errorfunction.hpp
@brief Error function, its complement and its inverse
*/

#ifndef QF_ERRORFUNCTION_HPP
#define QF_ERRORFUNCTION_HPP

#include <qflib/defines.hpp>
#include <qflib/exception.hpp>
#include <cmath>

BEGIN_NAMESPACE(qf)

/** Class containing functions related to the error function.
    The functions are static members, no objects of class ErrorFunction can be created.
    erf/erfc use the C++ standard library (std::erf, std::erfc).
    inverf/inverfc use a rational approximation with Halley refinement.
*/
class ErrorFunction
{
public:
  /** Returns the error function erf(x) */
  static double erf(double x);

  /** Returns the complement 1-erf(x) */
  static double erfc(double x);

  /** Returns the inverse of the error function */
  static double inverf(double p);

  /** Returns the inverse of the complement of the error function */
  static double inverfc(double p);

private:
  // This class is a namespace for error function related functions.
  // No default or copy construction or assignment allowed.
  ErrorFunction() = delete;
  ErrorFunction(ErrorFunction const&) = delete;
  ErrorFunction& operator=(ErrorFunction const&) = delete;
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline double ErrorFunction::erf(double x)
{
  return std::erf(x);
}

inline double ErrorFunction::erfc(double x)
{
  return std::erfc(x);
}

inline double ErrorFunction::inverf(double p)
{
  return inverfc(1. - p);
}

END_NAMESPACE(qf)

#endif // ORF_ERRORFUNCTION_HPP
