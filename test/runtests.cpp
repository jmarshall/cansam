#include <iostream>
#include <cstdlib>
#include <cstring>

#include "sam/exception.h"

#include "test/test.h"

void test_harness::check(bool expr, const string& title) {
  if (expr) {
    npass++;
  }
  else {
    std::cerr << "FAILED\t" << title << '\n';
    nfail++;
  }
}

void test_harness::check(const string& a, const string& b, const string& title){
  if (a == b) {
    npass++;
  }
  else {
    std::cerr << "FAILED\t" << title << " (\"" << a << "\" != \"" << b
	      << "\")\n";
    nfail++;
  }
}

void test_harness::check(const char* a, const char* b, const string& title) {
  if (strcmp(a, b) == 0) {
    npass++;
  }
  else {
    std::cerr << "FAILED\t" << title << " (\"" << a << "\" != \"" << b
	      << "\")\n";
    nfail++;
  }
}

void test_harness::check(size_t a, size_t b, const string& title) {
  if (a == b) {
    npass++;
  }
  else {
    std::cerr << "FAILED\t" << title << " (" << a << " != " << b << ")\n";
    nfail++;
  }
}

extern void test_alignments(test_harness&);
extern void test_headers(test_harness&);
extern void test_sam_io(test_harness&);
extern void test_wire(test_harness&);

int main() {
  test_harness t;

  try {
    test_headers(t);
    test_alignments(t);
    test_sam_io(t);
    test_wire(t);
  }
  catch (const sam::exception& e) {
    if (e.filename().empty())
      std::cerr << "Excess sam::exception: " << e.what() << '\n';
    else
      std::cerr << "Excess sam::exception from " << e.filename() << ": "
		<< e.what() << '\n';
    t.nfail++;
  }
  catch (const std::exception& e) {
    std::cerr << "Excess std::exception: " << e.what() << '\n';
    t.nfail++;
  }
  catch (const char* what) {
    std::cerr << "Excess exception: " << what << '\n';
    t.nfail++;
  }

  if (t.nfail > 0) {
    std::cerr << "\nTotal failures: " << t.nfail << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
