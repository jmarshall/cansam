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
  /// that of @a b.
  bool operator() (const alignment& a, const alignment& b) const
    { return cmp_by_qname(a, b) < 0; }
};

} // namespace sam

#endif
