#ifndef CANSAM_BGZF_H
#define CANSAM_BGZF_H

#include <cstring>
#include "lib/wire.h"

namespace sam {
namespace BGZF {

/* Utilities for identifying a BGZF block header, i.e., a GZIP member header
(RFC 1952) with a 'BC' extra subfield.  */

// Size, in bytes, of a BGZF block header
enum { hsize = 18 };

// Returns whether the specified memory block starts with a GZIP member header
inline bool is_gzip_header(const char* s, int length) {
  return length >= 2 && s[0] == '\x1f' && s[1] == '\x8b';
}

// Returns whether the specified memory block starts with a valid BGZF header
inline bool is_bgzf_header(const char* s, int length) {
  return length >= 18 && s[0] == '\x1f' && s[1] == '\x8b' &&
	 (s[3] & 4) && memcmp(&s[10], "\6\0\x42\x43\2\0", 6) == 0;
}

// For a valid BGZF header, returns the block_size field
inline int block_size(const char* text) {
  return convert::uint16(&text[16]) + 1;
}

} // namespace BGZF
} // namespace sam

#endif
