#ifndef CANSAM_UTILITIES_H
#define CANSAM_UTILITIES_H

#include <string>

#include "sam/types.h"

namespace sam {

// Removes a trailing line terminator, whether it be LF, CR, or CR-LF.
// (Usually CR would be because there was a CR-LF terminator and the LF has
// already been elided.)
std::string& chomp(std::string& s);

coord_t to_int(const std::string& str, std::string::size_type begin,
	       std::string::size_type end);

} // namespace sam

#endif
