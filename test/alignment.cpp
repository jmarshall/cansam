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

void test_iterators(test_harness& t, sam::alignment& aln) {
  sam::alignment::iterator it;
  it = aln.begin();
  sam::alignment::iterator it1 = it;
  sam::alignment::iterator it2(it);

  sam::alignment::const_iterator cit;
  cit = aln.begin();
  sam::alignment::const_iterator cit1 = cit;
  sam::alignment::const_iterator cit2(cit);

  t.check(it1 == it2, "iterator.==");
  t.check(!(it1 != it2), "iterator.!=");

  t.check(cit1 == cit2, "const_iterator.==");
  t.check(!(cit1 != cit2), "const_iterator.!=");

  t.check(it1 == cit2, "iterator__const_iterator.==");
  t.check(cit1 == it2, "const_iterator__iterator.==");
  t.check(!(it1 != cit2), "iterator__const_iterator.!=");
  t.check(!(cit1 != it2), "const_iterator__iterator.!=");

std::cout << "begin: " << aln.begin() << ", end: " << aln.end() << "\n";

  std::string foo = "foo";
  int bar = 37;

aln.dump_on(std::cout);
std::cout << aln << "\npush_back.X1\n";
  aln.push_back("X1", foo);
aln.dump_on(std::cout, aln.find("X1"));
std::cout << aln << "\npush_back.X2\n";
  aln.push_back("X2", bar);
  //aln.push_back("X3", cit2);

aln.dump_on(std::cout);
std::cout << aln << "\ninsert.Y1\n";
  aln.insert(aln.begin(), "Y1", foo);
aln.dump_on(std::cout, aln.begin());
std::cout << aln << "\ninsert.Y2\n";
  aln.insert(aln.begin(), "Y2", bar);
  //aln.insert(aln.begin(), "Y3", cit2);

  //aln.replace(aln.begin(), aln.end(), "Y3", it2);

aln.dump_on(std::cout);
std::cout << aln << "\nset_aux.Z1\n";
  aln.set_aux("Z1", foo);
aln.dump_on(std::cout);
std::cout << aln << "\nset_aux.Z2\n";
  aln.set_aux("Z2", bar);
aln.dump_on(std::cout);
std::cout << aln << "\nset_aux.Z3\n";
  aln.set_aux("Z3", cit2);

std::cout << "begin: " << aln.begin() << ", end: " << aln.end() << "\nauxen:";
  for (cit = aln.begin(); cit != aln.end(); ++cit)
    std::cout << " {" << *cit << "}";
std::cout << "\n";

  it1 = aln.begin(), ++it1, ++it1;
aln.dump_on(std::cout, it1);
std::cout << aln << "\nset_aux.it1\n";
  aln.set_aux(it1, foo);
aln.dump_on(std::cout, it1);
std::cout << aln << "\nset_aux.it1\n";
  aln.set_aux(it1, bar);
  //aln.set_aux(it1, cit2);

  it = aln.begin();
std::cout << "And finally...\n";
  aln.dump_on(std::cout, it);
std::cout << "pre:  it: " << it << "\n";
  it++;
std::cout << "post: it: " << it << "\n";
  aln.dump_on(std::cout, it);
std::cout << "End of test_iterators()\n";
}

void test_alignments(test_harness& t) {
  sam::alignment aln;

  t.check(sizeof aln, sizeof(void*), "aln.sizeof");

  sam::alignment a1(aln);
  sam::alignment a2 = aln;

  a2 = a1;

  a1.set_qname("JAS5_12:1:3");

  test_unpack_seq(t);
  test_iterators(t, a1);
}
