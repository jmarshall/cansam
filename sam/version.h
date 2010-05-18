/// @file sam/version.h
/// Cansam library version information

#ifndef CANSAM_VERSION_H
#define CANSAM_VERSION_H

#include <string>

/** @file
The Cansam library version number is a 3-part @c major.minor.patch number,
provided both as a preprocessor constant and a text string.  While the former
denotes the Cansam headers your code is compiled against and the latter the
Cansam library used at link- or run-time, in most circumstances both will
represent the same value.

The preprocessor macro #CANSAM_VERSION encodes the version as an integer
literal in the same way as <a href="http://www.boost.org/">Boost</a>'s
version macro, with two decimal digits of @c patch and three of @c minor.
The version() function returns the version number in text form, with the
final @c ".nn" omitted if @c patch is 0.

This header file is not included by other Cansam header files, so your code
should @c @#include it itself where necessary.

@note While the utilities supplied with the Cansam library use this library
version number directly as their own version number, external code should not
-- as such code is not updated in concert with the library.  */

/// Library version, as an integer @c Mmmmpp
#define CANSAM_VERSION  690

namespace sam {

/// Library version, as @c "n.nnn[.nn]"
std::string version();

} // namespace sam

#endif
