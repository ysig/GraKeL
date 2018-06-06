/* Arash Partov Function
 * Author: Ioannis Siglidis <y.siglidis@gmail.com>
 * License: BSD 3 clause"
 * Code taken from: http://www.partow.net/programming/hashfunctions/#APHashFunction
 */
#include "../include/functions.hpp"

unsigned int ArashPartov(const char* str, unsigned int length)
{
   unsigned int hash = 0xAAAAAAAA;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash ^= ((i & 1) == 0) ? (  (hash <<  7) ^ (*str) * (hash >> 3)) :
                               (~((hash << 11) + ((*str) ^ (hash >> 5))));
   }

   return hash;
}
