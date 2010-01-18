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

/// Returns an incremented copy of the given @a iterator
/** This is a simple version of the upcoming <tt>std::%next()</tt>
function template. */
template <typename InputIterator>
inline InputIterator next(InputIterator iterator) { return ++iterator; }

} // namespace sam

// @cond private
#define no_throw() throw()
// @endcond

#endif
