#include <iostream> // FIXME NUKE-ME

#include "test/test.h"
#include "sam/alignment.h"

std::string unpack_seq(const char* raw_seq, int seq_length) {
  std::string s;
  sam::alignment::unpack_seq(s, raw_seq, seq_length);
  return s;
}

void test_packing(test_harness& t,
		  string seq, string packed, string unpacked = string()) {
  char* buf = new char[packed.length()];

  if (unpacked.empty())  unpacked = seq;

  sam::alignment::pack_seq(buf, seq.c_str(), seq.length());
  t.check(string(buf, packed.length()), packed, "pack_seq." + seq);
  t.check(unpack_seq(buf, unpacked.length()), unpacked, "unpack_seq." + seq);

  delete [] buf;
}

void test_unpack_seq(test_harness& t) {
  t.check(unpack_seq("", 0), "", "unpack_seq.empty");
  t.check(unpack_seq("\x10", 1), "A", "unpack_seq.A");
  t.check(unpack_seq("\x20", 1), "C", "unpack_seq.C");
  t.check(unpack_seq("\x40", 1), "G", "unpack_seq.G");
  t.check(unpack_seq("\x80", 1), "T", "unpack_seq.T");
  t.check(unpack_seq("\x18", 2), "AT", "unpack_seq.AT");
  t.check(unpack_seq("\x84", 2), "TG", "unpack_seq.TG");
  t.check(unpack_seq("\x42", 2), "GC", "unpack_seq.GC");
  t.check(unpack_seq("\x41\x88\x12\x10", 7), "GATTACA", "unpack_seq.GATTACA");

  test_packing(t, "", "");
  test_packing(t, "AACGTCCTGCAGGATTA", "\x11\x24\x82\x28\x42\x14\x41\x88\x10");
  test_packing(t, "acmgrsvtwyhkdb.", "\x12\x34\x56\x78\x9a\xbc\xde\xf0",
		  "ACMGRSVTWYHKDBN");

  bool threw = false;
  try {
    char buffer[10];
    sam::alignment::pack_seq(buffer, "ATZ", 3);
  }
  catch (...) { threw = true; }
  t.check(threw, "pack_seq.invalid_char");
}

void test_iterators(test_harness& t, const sam::alignment& aln) {
  sam::alignment::const_iterator it;
  sam::alignment::const_iterator it1 = it;
  sam::alignment::const_iterator it2(it);

  t.check(it1 == it2, "const_iterator.==");
  t.check(!(it1 != it2), "const_iterator.!=");

std::clog << "begin: " << aln.begin() << ", end: " << aln.end() << "\n";

  for (it = aln.begin(); it != aln.end(); ++it) {
    // dereference it...
  }

  it = aln.begin();
std::clog << "And finally...\n";
  it++;
}

void test_alignments(test_harness& t) {
  sam::alignment aln;

  t.check(sizeof aln, sizeof(void*), "aln.sizeof");

  sam::alignment a1(aln);
  sam::alignment a2 = aln;

  a2 = a1;

  test_unpack_seq(t);
  test_iterators(t, a1);
}
