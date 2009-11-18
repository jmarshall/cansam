#include <iostream>

#include "test/test.h"
#include "lib/wire.h"

using namespace sam;

union wire16 {
  uint16_t u;
  char c[2];
  };

union wire32 {
  uint32_t u;
  char c[4];
  };

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
