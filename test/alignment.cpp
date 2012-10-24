/*  test/alignment.cpp -- Tests for alignment records.

    Copyright (C) 2010-2012 Genome Research Ltd.

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

#include <iostream> // FIXME NUKE-ME
#include <sstream>

#include "test/test.h"
#include "cansam/sam/alignment.h"

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
  sam::alignment aln2 = aln;
  cit2 = aln2.find("X1");
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

void test_auxen(test_harness& t) {
  sam::alignment aln;
  aln.push_back("XS", "carrot");
  aln.push_back("XI", 37);

  t.check(aln.aux<string>("XS"), "carrot", "aux<string>.1");
  t.check(aln.aux<string>("XI"), "37", "aux<string>.2");

  t.check(aln.aux<const char*>("XS"), "carrot", "aux<const char*>");
  t.check(aln.aux<int>("XI"), 37, "aux<int>");
}

void test_format(test_harness& t, std::ios::fmtflags fmt, const char* prefix) {
  char buffer[64];

  for (int i = 0; i <= 0x7ff; i++) {
    std::ostringstream s;
    s.setf(fmt, std::ios::basefield|std::ios::showbase|std::ios::uppercase);
    if ((fmt & std::ios::hex) && i != 0)  s << "0x";
    s << i;

    *sam::format_sam(buffer, i, s) = '\0';
    t.check(buffer, s.str(), prefix + s.str());
  }
}

void test_format(test_harness& t) {
  test_format(t, std::ios::dec, "dec.");
  test_format(t, std::ios::oct|std::ios::showbase,  "oct.");
  test_format(t, std::ios::hex|std::ios::uppercase, "hex.");

  char buffer[64];
  std::ostringstream s;
  s << std::boolalpha;

  for (int i = 0; i <= 0x7ff; i++) {
    *sam::format_sam(buffer, i, s) = '\0';
    t.check(sam::parse_flags(buffer), i, std::string("symbolic.") + buffer);
  }
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
  test_auxen(t);

  test_format(t);
}
