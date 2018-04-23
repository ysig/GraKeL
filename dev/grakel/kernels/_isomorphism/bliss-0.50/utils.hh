#ifndef BLISS_UTILS_HH
#define BLISS_UTILS_HH

/**
 * \file
 * \brief Some small utilities.
 *
 */

/*
 * Copyright (c) Tommi Junttila
 * Released under the GNU General Public License version 2.
 */

#include <cstdio>

namespace bliss {

/**
 * Print the permutation \a perm of {0,...,N-1} in the cycle format
 * in the file stream \a fp.
 * The amount \a offset is added to each element before printing,
 * e.g. the permutation (2 4) is printed as (3 5) when \a offset is 1.
 */
void print_permutation(FILE *fp,
		       const unsigned int N,
		       const unsigned int *perm,
		       const unsigned int offset = 0);

} // namespace bliss

#endif
