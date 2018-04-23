#ifndef BLISS_UINTSEQHASH_HH
#define BLISS_UINTSEQHASH_HH

#include <cstdio>

namespace bliss {

/**
 * \brief A hash for sequences of unsigned ints.
 */
class UintSeqHash
{
protected:
  unsigned int h;
public:
  UintSeqHash() {h = 0; }
  UintSeqHash(const UintSeqHash &other) {h = other.h; }
  UintSeqHash& operator=(const UintSeqHash &other) {h = other.h; return *this; }
  
  /** Reset the hash value. */
  void reset() {h = 0; }

  /** Add the unsigned int \a n to the sequence. */
  void update(unsigned int n);

  /** Get the hash value of the sequence seen so far. */
  unsigned int get_value() const {return h; }

  /** Compare the hash values of this and \a other.
   * Return -1/0/1 if the value of this is smaller/equal/greater than
   * that of \a other. */
  int cmp(const UintSeqHash &other) const {
    return (h < other.h)?-1:((h == other.h)?0:1);
  }
  /** An abbreviation for cmp(other) < 0 */
  bool is_lt(const UintSeqHash &other) const {return(cmp(other) < 0); }
  /** An abbreviation for cmp(other) <= 0 */
  bool is_le(const UintSeqHash &other) const {return(cmp(other) <= 0); }
  /** An abbreviation for cmp(other) == 0 */
  bool is_equal(const UintSeqHash &other) const {return(cmp(other) == 0); }
};


} // namespace bliss

#endif
