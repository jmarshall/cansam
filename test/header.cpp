#include "sam/header.h"
#include "test/test.h"

using sam::header;

void test_headers(test_harness& t) {
  const char sq[] = "@SQ	SN:foo	LN:15	SP:human";
  const char pg[] = "@PG	ID:runtests	VN:3.14	CL:test/runtests";

  header h(sq);
  t.check(h.str(), sq, "sq: ctor/str");

  h.assign(pg);
  t.check(h.str(), pg, "pg: assign/str");

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
}
