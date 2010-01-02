#ifndef CANSAM_UTILITIES_H
#define CANSAM_UTILITIES_H

#include <string>
#include <sstream>

#include "sam/types.h"

namespace sam {

// Removes a trailing line terminator, whether it be LF, CR, or CR-LF.
// (Usually CR would be because there was a CR-LF terminator and the LF has
// already been elided.)
std::string& chomp(std::string& s);

coord_t to_int(const std::string& str, std::string::size_type begin,
	       std::string::size_type end);

class make_string {
public:
  make_string() { }

  make_string& operator<< (const char* text) { buffer << text; return *this; }
  make_string& operator<< (char c) { unsigned char uc = c; return *this << uc; }
  make_string& operator<< (unsigned char c);

#if 0
  template<typename T>
  make_string& operator<< (const T& t) { buffer << t; return *this; }
#endif

  operator std::string () const { return buffer.str(); }

private:
  std::ostringstream buffer;
};

} // namespace sam

#endif
