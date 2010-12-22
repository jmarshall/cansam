/*  test/wire.cpp -- Tests for binary data access routines.

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

#include <iostream>

#include "test/test.h"
#include "lib/wire.h"

using namespace sam::convert;

union wire16 {
  uint16_t u;
  char c[2];
  };

union wire32 {
  uint32_t u;
  char c[4];
  };

#if 0
// FIXME Test *after* lib/wire.h settles down...

void test_wire(test_harness& t) {
  std::ios::fmtflags flags = std::cerr.flags();
  std::cerr << std::hex << std::showbase;

  wire16 u16;
  u16.c[1] = 0x12, u16.c[0] = 0x34;
  t.check(bamtoh16(u16.u), 0x1234, "bamtoh16.1234");

  u16.c[1] = 0xAB, u16.c[0] = 0xCD;
  t.check(bamtoh16(u16.u), 0xABCD, "bamtoh16.ABCD");

  wire32 u32;
  u32.c[3] = 0x12, u32.c[2] = 0x34, u32.c[1] = 0x56, u32.c[0] = 0x78;
  t.check(bamtoh32(u32.u), 0x12345678, "bamtoh32.12345678");

  u32.c[3] = 0x89, u32.c[2] = 0xAB, u32.c[1] = 0xCD, u32.c[0] = 0xEF;
  t.check(bamtoh32(u32.u), 0x89ABCDEF, "bamtoh32.89ABCDEF");

  t.check(bamtoh16(htobam16(0x89EF)), 0x89EF, "idem16.89EF");
  t.check(bamtoh16(htobam16(0x4321)), 0x4321, "idem16.4321");
  t.check(bamtoh32(htobam32(0x89ABCDEF)), 0x89ABCDEF, "idem32.89ABCDEF");
  t.check(bamtoh32(htobam32(0x76543210)), 0x76543210, "idem32.76543210");

  std::cerr.flags(flags);
}

#else
void test_wire(test_harness&) { }
#endif
