#include "test/test.h"
#include "sam/alignment.h"

static void test_unpack_seq(test_harness& t) {
  t.check(sam::alignment::unpack_seq("", 0), "", "unpack_seq.empty");
  t.check(sam::alignment::unpack_seq("\x10", 1), "A", "unpack_seq.A");
  t.check(sam::alignment::unpack_seq("\x20", 1), "C", "unpack_seq.C");
  t.check(sam::alignment::unpack_seq("\x40", 1), "G", "unpack_seq.G");
  t.check(sam::alignment::unpack_seq("\x80", 1), "T", "unpack_seq.T");
  t.check(sam::alignment::unpack_seq("\x18", 2), "AT", "unpack_seq.AT");
  t.check(sam::alignment::unpack_seq("\x84", 2), "TG", "unpack_seq.TG");
  t.check(sam::alignment::unpack_seq("\x42", 2), "GC", "unpack_seq.GC");
  t.check(sam::alignment::unpack_seq("\x41\x88\x12\x10", 7), "GATTACA",
	  "unpack_seq.GATTACA");
}

void test_alignments(test_harness& t) {
  sam::alignment aln;

  t.check(sizeof aln, sizeof(void*), "aln.sizeof");

  test_unpack_seq(t);
}
