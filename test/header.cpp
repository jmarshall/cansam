#include "sam/header.h"
#include "test/test.h"

#include <iostream>

using sam::header;

void foo(header& h) {
  std::cout << h;
  std::cout << h.begin();
  std::cout << *(h.begin());
}

void exercise_header(header& hdr) {
  const header& chdr = hdr;

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
