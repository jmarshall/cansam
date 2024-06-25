/*  wire.h -- Access binary data irrespective of endianness and alignment.

    Copyright (C) 2010 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
 3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
    nor the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND ITS CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */

#ifndef WIRE_H
#define WIRE_H

#include <stdint.h>

// Integers in BAM files are stored in a little-endian format.  We use the
// macros predefined by the various compilers to detect what needs to be done
// to convert such a binary representation to that used by this host.

#if defined __i386__ || defined __x86_64__ || defined __ia64__ || \
    defined __arm__ || defined __aarch64__
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

/* Conversion functions are provided for various combinations of aligned vs.
unaligned memory, and converting in-place vs. reading-as-host/writing-as-bam.

Convert a variable or aligned memory in-place...
   void set_uint16(uint16_t&)     ...to host representation
   void set_bam16(uint16_t&)      ...to the BAM wire format

Convert unaligned memory in-place...
   void set_bam16(void*)          ...to the BAM wire format

Return a (host-representation) value...
   uint16_t uint16(const void*)   ...by reading BAM-formatted unaligned memory

Write a (host-represented) value...
   void set_bam_uint16(void*, uint16_t)  ...to unaligned memory in BAM format

Similar functions are provided for {int,uint}{16,32}.  */
// FIXME Also likely eventually for {int,uint}{64}.

#if defined WIRE_NOOP

inline void set_uint16(uint16_t&) { }
inline void set_uint32(uint32_t&) { }
inline void set_bam16(uint16_t&) { }
inline void set_bam32(uint32_t&) { }
inline void set_bam16(void*) { }
inline void set_bam32(void*) { }

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

inline void set_int16(int16_t& x) { uint16_t ux = x; set_uint16(ux); x = ux; }
inline void set_int32(int32_t& x) { uint32_t ux = x; set_uint32(ux); x = ux; }
inline void set_bam16(int16_t& x) { uint16_t ux = x; set_bam16(ux); x = ux; }
inline void set_bam32(int32_t& x) { uint32_t ux = x; set_bam32(ux); x = ux; }

inline int16_t int16(const void* pv) { return uint16(pv); }
inline int32_t int32(const void* pv) { return uint32(pv); }
inline void set_bam_int16(void* pv, int16_t x) { set_bam_uint16(pv, x); }
inline void set_bam_int32(void* pv, int32_t x) { set_bam_uint32(pv, x); }

} // namespace convert
} // namespace sam

#endif
