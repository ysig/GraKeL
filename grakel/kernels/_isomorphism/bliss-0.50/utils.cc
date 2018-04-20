#include <assert.h>
#include "utils.hh"

namespace bliss {

void print_permutation(FILE *fp,
		       const unsigned int N,
		       const unsigned int *perm,
		       const unsigned int offset)
{
  assert(N > 0);
  assert(perm);
  for(unsigned int i = 0; i < N; i++) {
    unsigned int j = perm[i];
    if(j == i)
      continue;
    bool is_first = true;
    while(j != i) {
      if(j < i) {
        is_first = false;
        break;
      }
      j = perm[j];
    }
    if(!is_first)
      continue;
    fprintf(fp, "(%u,", i+offset);
    j = perm[i];
    while(j != i) {
      fprintf(fp, "%u", j+offset);
      j = perm[j];
      if(j != i)
        fprintf(fp, ",");
    }
    fprintf(fp, ")");
  }
}

} // namespace bliss
