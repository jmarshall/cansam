#ifndef WIRE_H
#define WIRE_H

#include <stdint.h>

// Integers in BAM files are stored in a little-endian format.  We use the
// macros predefined by the various compilers to detect what needs to be done
// to convert such a binary representation to that used by this host.

#if defined __i386__ || defined __x86_64__ || defined __ia64__
#define WIRE_NOOP
#elif 0
// We need to swap bytes, and the architecture does rotates efficiently.
#define WIRE_ROTATE_BYTES
#else
#error - Unknown architecture
#endif

/* Converting integers to and from a wire format such as the BAM format
requires solutions to two problems: translation between the host and wire
representations, which usually means converting endianness; and, unless
the format is particularly well-designed, accessing wire-format items that
are not properly aligned.
*/

namespace sam {
namespace convert {

/* FIXME
Possible signatures:
  aligned memory (so value, not ptr)  v  unaligned memory (ptr)
  swap in place  v  read-as-host/write-as-bam

Swap in place:  (probably only unaligned (is more general) is actually needed)
  void swap(uint16_t*)	-- un/hostify aligned memory in place
  void swap(uint16_t&)	-- un/hostify register or aligned memory in place
? void swap(void*)	-- un/hostify unaligned memory in place

(actually unaligned swap is useless as it doesn't solve the how-will-you-look-
at-it---it's-unaligned problem)

Read as host:
? uint16_t hostify(uint16_t)     -- aligned memory, no explicit ptr needed
# uint16_t hostify(const void*)  -- unaligned memory

Write as bam:
? uint16_t bamify(uint16_t) -- aligned memory, no explicit ptr needed
  void bamify(void* dest, uint16_t) -- aligned memory (probably not useful)
# void bamify(void* dest, uint16_t) -- unaligned memory

Proposed names:

? void set_uint16(uint16_t&)	-- convert unaligned memory bam to host
? void set_bam16(uint16_t&) 	-- convert unaligned memory host to bam

? uint16_t uint16(uint16_t)	-- convert bam to host (for aligned)
# uint16_t uint16(const void*)	-- read unaligned memory as host

? uint16_t bam_uint16(uint16_t) 	-- convert host to bam (for aligned)
# void set_bam_uint16(void*, uint16_t)	-- write bam to unaligned memory

# = already in use
? = maybe useful to have this one
*/

#if defined WIRE_NOOP

#if 0
// FIXME Later, if they're needed...
inline void set_uint16(uint16_t&) { }
inline void set_uint32(uint32_t&) { }
inline void set_bam16(uint16_t&) { }
inline void set_bam32(uint32_t&) { }

inline uint16_t uint16(uint16_t x) { return x; }
inline uint32_t uint32(uint32_t x) { return x; }
inline uint16_t bam_uint16(uint16_t x) { return x; }
inline uint32_t bam_uint32(uint32_t x) { return x; }

inline int16_t int16(int16_t x) { return x; }
inline int32_t int32(int32_t x) { return x; }
inline int16_t bam_int16(int16_t x) { return x; }
inline int32_t bam_int32(int32_t x) { return x; }
#endif

// These are endianness-independent, and make the most pessimistic assumptions
// about not being able to do better than char-access to unaligned memory.
// (TODO On Intel, I think it's only bad if it spans the end of a cache line,
// and it may be better to do word-access 90% of the time even at the cost
// of a test.)

inline uint16_t uint16(const void* pv) {
  const unsigned char* p = static_cast<const unsigned char*>(pv);
  uint16_t x0 = p[0], x1 = p[1];
  return (x1 << 8) | x0;
}

inline uint32_t uint32(const void* pv) {
  const unsigned char* p = static_cast<const unsigned char*>(pv);
  uint32_t x0 = p[0], x1 = p[1], x2 = p[2], x3 = p[3];
  x3 = (x3 << 8) | x2;
  x3 = (x3 << 8) | x1;
  x3 = (x3 << 8) | x0;
  return x3;
}

inline void set_bam_uint16(void* pv, uint16_t x) {
  unsigned char* p = static_cast<unsigned char*>(pv);
  p[0] = x, x >>= 8;
  p[1] = x;
}

inline void set_bam_uint32(void* pv, uint32_t x) {
  unsigned char* p = static_cast<unsigned char*>(pv);
  p[0] = x, x >>= 8;
  p[1] = x, x >>= 8;
  p[2] = x, x >>= 8;
  p[3] = x;
}

inline int16_t int16(const void* pv) { return uint16(pv); }
inline int32_t int32(const void* pv) { return uint32(pv); }
inline void set_bam_int16(void* pv, int16_t x) { set_bam_uint16(pv, x); }
inline void set_bam_int32(void* pv, int32_t x) { set_bam_uint32(pv, x); }

#elif defined WIRE_ROTATE_BYTES

inline void set_uint16(void* pv) {
  char* p = static_cast<char*>(pv);
  char tmp;
  tmp = p[0], p[0] = p[1], p[1] = tmp;
}

inline void set_uint32(void* pv) {
  char* p = static_cast<char*>(pv);
  char tmp;
  tmp = p[0], p[0] = p[3], p[3] = tmp;
  tmp = p[1], p[1] = p[2], p[2] = tmp;
}

inline void set_bam16(void* pv) { set_uint16(pv); }
inline void set_bam32(void* pv) { set_uint32(pv); }

inline uint16_t uint16(uint16_t x) { return (x << 8) | (x >> 8); }

inline uint32_t uint32(uint32_t x)
  { uint32_t xA_C_ = x & 0xff00ff00;
    uint32_t x_B_D = x & 0x00ff00ff;
    uint32_t x_C_A = (xA_C_ << 8) | (xA_C_ >> 24);
    uint32_t xD_B_ = (x_B_D >> 8) | (x_B_D << 24);
    return xD_B_ | x_C_A; }

inline uint16_t bam_uint16(uint16_t x) { return uint16(x); }
inline uint32_t bam_uint32(uint32_t x) { return uint32(x); }

#endif

} // namespace convert
} // namespace sam

#endif
