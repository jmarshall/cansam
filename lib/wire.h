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
#error unknown architecture
#endif

namespace sam {

/* FIXME
Possible signatures:
  aligned memory (so value, not ptr)  v  unaligned memory (ptr)
  swap in place  v  read-as-host/write-as-bam

Swap in place:  (probably only unaligned (is more general) is actually needed)
  void swap(uint16_t*)	-- un/hostify aligned memory in place
  void swap(void*)	-- un/hostify unaligned memory in place

Read as host:
  uint16_t hostify(uint16_t)     -- aligned memory, no explicit ptr needed
  uint16_t hostify(const void*)  -- unaligned memory

Write as bam:
  uint16_t bamify(uint16_t) -- aligned memory, no explicit ptr needed
  void bamify(void* dest, uint16_t) -- aligned memory (probably not useful)
  void bamify(void* dest, uint16_t) -- unaligned memory
*/

#if defined WIRE_NOOP

inline void hostify16(void*) { }
inline void bamify16(void*) { }

inline void cvt_bamtoh32(void*) { }
inline void cvt_htobam32(void*) { }

inline uint16_t bamtoh16(uint16_t x) { return x; }
inline uint16_t htobam16(uint16_t x) { return x; }

inline uint32_t bamtoh32(uint32_t x) { return x; }
inline uint32_t htobam32(uint32_t x) { return x; }

#elif defined WIRE_ROTATE_BYTES

inline void cvt_bamtoh16(void* pv)
  { char* p = static_cast<char*>(pv);
    char tmp = p[0]; p[0] = p[1]; p[1] = tmp; }

inline void cvt_htobam16(void* pv) { cvt_bamtoh16(pv); }

inline void cvt_bamtoh32(void* pv)
  { char* p = static_cast<char*>(pv);
    char tmp;
    tmp = p[0]; p[0] = p[3]; p[3] = tmp;
    tmp = p[1]; p[1] = p[2]; p[2] = tmp; }

inline void cvt_htobam32(void* pv) { cvt_bamtoh32(pv); }

inline uint16_t bamtoh16(uint16_t x) { return (x << 8) | (x >> 8); }
inline uint16_t htobam16(uint16_t x) { return (x << 8) | (x >> 8); }

inline uint32_t bamtoh32(uint32_t x)
  { uint32_t xA_C_ = x & 0xff00ff00;
    uint32_t x_B_D = x & 0x00ff00ff;
    uint32_t x_C_A = (xA_C_ << 8) | (xA_C_ >> 24);
    uint32_t xD_B_ = (x_B_D >> 8) | (x_B_D << 24);
    return xD_B_ | x_C_A; }

inline uint32_t htobam32(uint32_t x) { return bamtoh32(x); }

#endif

} // namespace sam

#endif
