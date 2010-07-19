/** @file   sam/algorithm.h
    @brief  Algorithms operating on alignment records
*/

#ifndef CANSAM_ALGORITHM_H
#define CANSAM_ALGORITHM_H

#include <functional>

#include "sam/alignment.h"

namespace sam {

/// Function object class for comparing alignments
class less_by_qname :
  public std::binary_function<const alignment&, const alignment&, bool> {
public:
  /// Returns whether the query name of @a a is lexicographically less than
  /// that of @a b
  bool operator() (const alignment& a, const alignment& b) const
    { return cmp_by_qname(a, b) < 0; }
};

/// Function object class for comparing alignments
class equal_by_qname :
  public std::binary_function<const alignment&, const alignment&, bool> {
public:
  /// Returns whether the query name of @a a is the same as that of @a b
  bool operator() (const alignment& a, const alignment& b) const
    { return cmp_by_qname(a, b) == 0; }
};

/// Function object class for hashing alignments
class hash_by_qname : public std::unary_function<const alignment&, size_t> {
public:
  /// Returns a hash value derived from the alignment's query name
  size_t operator() (const alignment& aln) const {
    size_t result = 0;
    for (const char* s = aln.qname_c_str(); *s; s++)
      result = (result * 131) + *s;
    return result;
  }
};

} // namespace sam

#endif
