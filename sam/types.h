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

} // namespace sam

#endif
