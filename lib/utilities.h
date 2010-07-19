#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <sstream>
#include <limits>

#include "sam/types.h"
#include "lib/bits/sign_traits.h"

namespace sam {

/* These functions format the given VALUE into the DEST buffer and return
a pointer to the first unused buffer position.  The constants give the
minimum buffer size that should be used with the corresponding functions.  */
// FIXME doco
namespace format {

template <typename IntType>
struct buffer {
  /* digits10 is the "number of base 10 digits that can be represented without
  change", which is usually one less than the maximum number of base 10 digits
  that can be represented -- e.g., consider 9999 v 65535.  Thus we add 2, one
  for the excess digit and one for a sign character.  */
  static const int size = std::numeric_limits<IntType>::digits10 + 2;

  // TODO Add buffer sizes for octal and hexadecimal
};

template <typename UnsignedType>
char* decimal_(char* dest, UnsignedType value, const traits::false_type&) {
  UnsignedType n = value;
  do { dest++; n /= 10; } while (n != 0);

  char* destlim = dest;
  do { *--dest = (value % 10) + '0'; value /= 10; } while (value != 0);

  return destlim;
}

template <typename SignedType>
char* decimal_(char* dest, SignedType value, const traits::true_type&) {
  typename traits::make_unsigned<SignedType>::type uvalue = value;
  if (value < 0)  *dest++ = '-', uvalue = -uvalue;
  return decimal(dest, uvalue);
}

template <typename IntType>
char* decimal(char* dest, IntType value) {
  return decimal_(dest, value, traits::is_signed<IntType>());
}

template <typename IntType>
char* octal(char* dest, IntType ivalue) {
  typedef typename traits::make_unsigned<IntType>::type UnsignedType;
  UnsignedType value = ivalue;

  *dest++ = '0';

  UnsignedType n = value;
  while (n != 0) { dest++; n >>= 3; }

  char* destlim = dest;
  while (value != 0) { *--dest = (value & 07) + '0'; value >>= 3; }

  return destlim;
}

extern const char hexadecimal_digits[];

template <typename IntType>
char* hexadecimal(char* dest, IntType ivalue) {
  typedef typename traits::make_unsigned<IntType>::type UnsignedType;
  UnsignedType value = ivalue;

  *dest++ = '0';
  if (value != 0)  *dest++ = 'x';

  UnsignedType n = value;
  while (n != 0) { dest++; n >>= 4; }

  char* destlim = dest;
  while (value != 0) { *--dest = hexadecimal_digits[value & 0xf]; value >>= 4; }

  return destlim;
}

} // namespace format

namespace parse {

template <typename UnsignedType>
const char*
decimal_(const char* s, UnsignedType& value, const traits::false_type&) {
  value = 0;
  if (*s == '+')  s++;
  while (*s >= '0' && *s <= '9')  value = 10 * value + *s++ - '0';
  return s;
}

template <typename SignedType>
const char*
decimal_(const char* s, SignedType& value, const traits::true_type&) {
  typename traits::make_unsigned<SignedType>::type uvalue;
  if (*s == '-')  s = decimal(s+1, uvalue), value = -uvalue;
  else  s = decimal(s, uvalue), value = uvalue;
  return s;
}

template <typename IntType>
const char* decimal(const char* s, IntType& value) {
  return decimal_(s, value, traits::is_signed<IntType>());
}

} // namespace parse

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
