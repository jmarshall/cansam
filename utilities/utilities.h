#ifndef UTILITIES_UTILITIES_H
#define UTILITIES_UTILITIES_H

#include <iosfwd>

// Prints Cansam version number as PROGNAME's version number and
// brief copyright and (lack of) warranty information to STREAM.
void print_version(std::ostream& stream, const char* progname);

#endif
