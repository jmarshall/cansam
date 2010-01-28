#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <sstream>
#include <limits>

#if 0
#include <type_traits>
#else
#include <boost/type_traits/is_signed.hpp>
#include <boost/type_traits/make_unsigned.hpp>
#endif

#include "sam/types.h"

namespace sam {

/* These functions format the given VALUE into the DEST buffer and return
a pointer to the first unused buffer position.  The constants give the
minimum buffer size that should be used with the corresponding functions.
(digits10 is the "number of base 10 digits that can be represented without
change", which is usually one less than the maximum number of base 10 digits
that can be represented -- e.g., consider 9999 v 65535.  Thus we add 2, one
for the excess digit and one for a sign character.)  */
namespace format {

static const int int_digits = std::numeric_limits<int>::digits10 + 2;

template <typename UnsignedType>
char* decimal_(char* dest, UnsignedType value, const boost::false_type&) {
  UnsignedType n = value;
  do { dest++; n /= 10; } while (n != 0);

  char* destlim = dest;
  do { *--dest = (value % 10) + '0'; value /= 10; } while (value != 0);

  return destlim;
}

template <typename SignedType>
char* decimal_(char* dest, SignedType value, const boost::true_type&) {
  typename boost::make_unsigned<SignedType>::type uvalue = value;
  if (value < 0)  *dest++ = '-', uvalue = -uvalue;
  return decimal(dest, uvalue);
}

template <typename IntType>
char* decimal(char* dest, IntType value) {
  return decimal_(dest, value, boost::is_signed<IntType>());
}

} // namespace format

// Removes a trailing line terminator, whether it be LF, CR, or CR-LF.
// (Usually CR would be because there was a CR-LF terminator and the LF has
// already been elided.)
std::string& chomp(std::string& s);

coord_t to_int(const std::string& str, std::string::size_type begin,
	       std::string::size_type end);

class make_string {
public:
  make_string() { }

  template <typename T>
  make_string& operator<< (const T& t) { buffer << t; return *this; }

  make_string& operator<< (char c) { unsigned char uc = c; return *this << uc; }
  make_string& operator<< (unsigned char c);

  operator std::string () const { return buffer.str(); }

private:
  std::ostringstream buffer;
};

} // namespace sam

#endif
