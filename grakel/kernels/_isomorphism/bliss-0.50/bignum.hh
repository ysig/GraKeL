#ifndef BLISS_BIGNUM_HH
#define BLISS_BIGNUM_HH

/*
 * Copyright (c) Tommi Junttila
 * Released under the GNU General Public License version 2.
 */

#if defined(BLISS_USE_GMP)
#include <gmp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include "defs.hh"

namespace bliss {

/**
 * \brief A very simple class for big integers (or approximation of them).
 *
 * If the compile time flag BLISS_USE_GMP is set,
 * then the GNU Multiple Precision Arithmetic library (GMP) is used to
 * obtain arbitrary precision, otherwise "long double" is used to
 * approximate big integers.
 */


#if defined(BLISS_USE_GMP)


class BigNum
{
  mpz_t v;
public:
  /**
   * Create a new big number and set it to zero.
   */
  BigNum() {mpz_init(v); }

  /**
   * Destroy the number.
   */
  ~BigNum() {mpz_clear(v); }

  /**
   * Set the number to 'n'.
   */
  void assign(const int n) {mpz_set_si(v, n); }

  /**
   * Multiply the number with 'n'.
   */
  void multiply(const int n) {mpz_mul_si(v, v, n); }

  /**
   * Print the number in the file stream 'fp'.
   */
  int print(FILE *fp) {return mpz_out_str(fp, 10, v); }
};

#else

class BigNum
{
  long double v;
public:
  /**
   * Create a new big number and set it to zero.
   */
  BigNum(): v(0.0) {}

  /**
   * Set the number to 'n'.
   */
  void assign(const int n) {v = (long double)n; }

  /**
   * Multiply the number with 'n'.
   */
  void multiply(const int n) {v *= (long double)n; }

  /**
   * Print the number in the file stream 'fp'.
   */
  int print(FILE *fp) {return fprintf(fp, "%Lg", v); }
};

#endif

} //namespace bliss

#endif
