#include <iostream>
#include <sstream>

#include "sam/stream.h"
#include "test/test.h"

static void test_reader(test_harness& t) {
  std::istringstream sam1(
"foo\t37\t*\t0\t*\t*\t*\t0\t0\tATGC\t????\tNM:i:4\nbar\n");

  sam::isamstream str(sam1.rdbuf());
  sam::alignment aln;
  while (str >> aln)
    std::cout << aln << '\n';
}

void test_sam_io(test_harness& t) {
  test_reader(t);
}
