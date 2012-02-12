/*  test/header.cpp -- Tests for header records.

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

#include "cansam/sam/header.h"
#include "test/test.h"

#include <iostream>

using std::string;
using sam::header;

void test_foo(header& h) {
  std::cout << h;
  std::cout << h.begin();
  std::cout << *(h.begin());
}

void exercise_header(header& hdr) {
  const header& chdr = hdr;

  hdr.set_field("X1", 37);
  hdr.field<int>("X1");
  hdr.field("X1", 1);

  header::iterator it = hdr.begin();
  header::iterator& rit1 = ++it;
  header::iterator& rit2 = --it;
  const header::iterator& rit3 = it++;
  const header::iterator& rit4 = it--;
  header::iterator it1;
  it1 = it;
  header::iterator it2 = it;

  rit1 == rit2;
  rit3 != rit4;
  
  hdr.end();
  hdr.find("AA");
  hdr.empty();
  hdr.push_back("AB", "hello");
  hdr.push_back("AC", 37);
  hdr.insert(hdr.begin(), "AD", "hello");
  hdr.insert(hdr.begin(), "AE", 37);
  hdr.erase(hdr.begin());
  hdr.erase(hdr.begin(), hdr.end());
  hdr.clear();
  hdr.replace(hdr.begin(), hdr.end(), "AF", "hello");
  hdr.replace(hdr.begin(), hdr.end(), "AG", 37);

  header::const_iterator cit = chdr.begin();
  ++cit, --cit;
  cit++, cit--;
  header::const_iterator cit1;
  cit1 = cit;
  header::const_iterator cit2 = it;

  cit1 == cit2;
  cit1 != cit2;

  chdr.end();
  chdr.find("AA");
  chdr.empty();

  it == cit;
  it != cit;

  cit == it;
  cit != it;

  string s = cit->value<string>();
  int i = cit->value<int>();
  sam::coord_t x = cit->value<sam::coord_t>();
}

void test_headers(test_harness& t) {
#if 0
  const char sq[] = "@SQ	SN:foo	LN:15	SP:human";
  const char pg[] = "@PG	ID:runtests	VN:3.14	CL:test/runtests";

  header h(sq);
  t.check(h.str(), sq, "sq: ctor/str");

  h.assign(pg);
  t.check(h.str(), pg, "pg: assign/str");
#endif
  t.check(true, "happy");
#if 0

  h.assign("@NL");
  t.check(h.type(), "NL", "empty.type");
  t.check(h.size(), 0, "empty.size");

  h.assign("@CO	foobar");
  t.check(h.type(), "CO", "one.type");
  t.check(h.size(), 1, "one.size");
  t.check(h[0].value, "foobar", "one.value");

  h.assign(pg);
  t.check(h.size(), 3, "pg.size");
  t.check(h[0].tag, "ID", "pg[0].tag");
  t.check(h[0].value, "runtests", "pg[0].value");
  t.check(h[1].tag, "VN", "pg[1].tag");
  t.check(h[1].value, "3.14", "pg[1].value");
  t.check(h[2].tag, "CL", "pg[1].tag");
  t.check(h[2].value, "test/runtests", "pg[2].value");

  h.push_back("EX", "banana");
  t.check(h.size(), 4, "banana.size");
  t.check(h[3].tag, "EX", "banana[3].tag");
  t.check(h[3].value, "banana", "banana[3].value");

  t.check(h.field("ID"), "runtests", "field.found");
  t.check(h.field("XX", "missing"), "missing", "field.notfound");

  t.check(h.find("XX") == h.end(), "find.notfound");
#endif
}
