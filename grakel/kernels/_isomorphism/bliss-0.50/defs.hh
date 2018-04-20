#ifndef BLISS_DEFS_HH
#define BLISS_DEFS_HH

#include <cassert>

namespace bliss {

/**
 * The version number of bliss.
 */
static const char * const version = "0.50";


#if defined(BLISS_DEBUG)
#define BLISS_CONSISTENCY_CHECKS
#define BLISS_EXPENSIVE_CONSISTENCY_CHECKS
#endif

#if defined(BLISS_CONSISTENCY_CHECKS)
#define BLISS_ASSERT(a) assert(a)
//inline void BLISS_ASSERT(const int c) {assert(c); }
#else
#define BLISS_ASSERT(a) ;
//inline void BLISS_ASSERT(const int c) {}
#endif


#if defined(BLISS_CONSISTENCY_CHECKS)
/* Force a check that the found automorphisms are valid */
#define BLISS_VERIFY_AUTOMORPHISMS
#endif


#if defined(BLISS_CONSISTENCY_CHECKS)
/* Force a check that the generated partitions are equitable */
#define BLISS_VERIFY_EQUITABLEDNESS
#endif


} // namespace bliss



/*! \mainpage Bliss
 *
 * \section intro_sec Introduction
 *
 * This is the source code documentation of bliss,
 * produced by running <A href="http://www.doxygen.org">doxygen</A> in
 * the source directory.
 * The algorithms and data structures used in bliss are documented in
 * the papers found at the
 * <A href="http://www.tcs.hut.fi/Software/bliss">bliss web site</A>.
 *
 *
 * \section compile_sec Compiling
 *
 * Compiling bliss in Linux should be easy, just execute
 * \code
 * make
 * \endcode
 * in the bliss source directory.
 * This will produce the executable program \c bliss as well as
 * the library file \c libbliss.a that can be linked in other programs.
 * If you have the <A href="http://gmplib.org/">GNU Multiple Precision
 * Arithmetic Library</A> (GMP) installed in your machine, you can also use
 * \code
 * make gmp
 * \endcode
 * to enable exact computation of automorphism group sizes.
 *
 * When linking the bliss library \c libbliss.a in other programs,
 * remember to include the standard c++ library
 * (and the GMP library if you compiled bliss to include it).
 * For instance,
 * \code gcc -o test test.c -lstdc++ -lgmp -lbliss\endcode
 *
 * \section cppapi_sec The C++ language API
 *
 * The C++ language API is the main API to bliss;
 * all other APIs are just more or less complete variants of it.
 * The C++ API consists basically of the public methods in
 * the classes bliss::AbstractGraph, bliss::Graph, and bliss::Digraph.
 * For an example of its use,
 * see the \ref executable "source of the bliss executable".
 *
 *
 * \section capi_sec The C language API
 *
 * The C language API is given in the file bliss_C.h.
 * It is currently more restricted than the C++ API so
 * consider using the C++ API whenever possible.
 */


#endif
