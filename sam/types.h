/** @file  sam/types.h
    @brief Basic types for some SAM/BAM fields
*/

#ifndef CANSAM_TYPES_H
#define CANSAM_TYPES_H

namespace sam {

/// Type for positions and coordinates
typedef long coord_t;

/// @b Signed type for insert sizes, etc, or generally a difference
/// between coordinates
typedef long scoord_t;

#if 0
// @cond infrastructure
template <typename ValueType>
class indirect_iterator {
public:
  indirect_iterator& operator++ () { ++pit; return *this; }
  indirect_iterator operator++ (int)
    { indirect_iterator orig = *this; ++pit; return orig; }

  indirect_iterator& operator-- () { --pit; return *this; }
  indirect_iterator operator-- (int)
    { indirect_iterator orig = *this; --pit; return orig; }

private:
  UnderlyingIterator pit;
};
// @endcond
#endif

} // namespace sam

#endif
