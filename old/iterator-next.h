/** @file  sam/types.h
    @brief Basic types for some SAM/BAM fields
*/

#ifndef CANSAM_TYPES_H
#define CANSAM_TYPES_H

namespace sam {

/// Returns an incremented copy of the given @a iterator
/** This is a simple version of the upcoming <tt>std::%next()</tt>
function template.  */
template <typename InputIterator>
inline InputIterator next(InputIterator iterator) { return ++iterator; }

} // namespace sam

#endif
