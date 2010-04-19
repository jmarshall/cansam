#ifndef TEST_H
#define TEST_H

#include <string>
using std::string;

class test_harness {
public:
  test_harness() : npass(0), nfail(0) { }

  void check(bool expr, const string& title);
  void check(const string& a, const string& b, const string& title);
  void check(const char* a, const char* b, const string& title);
  void check(size_t a, size_t b, const string& title);

  int npass, nfail;
};

extern string test_objdir_prefix, test_srcdir_prefix;

#endif
