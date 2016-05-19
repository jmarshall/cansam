/*  test/runtests.cpp -- Simple test harness.

    Copyright (C) 2010-2012, 2016 Genome Research Ltd.

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

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "cansam/exception.h"

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

void test_harness::check(const std::ios& s, const string& title) {
  if (s) {
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
extern void test_intervals(test_harness&);
extern void test_sam_io(test_harness&);
extern void test_wire(test_harness&);

string test_objdir_prefix;
string test_srcdir_prefix;

int main(int argc, char** argv) {
  test_harness t;

  string slash = "/";
  test_objdir_prefix = (argc >= 2)? argv[1] + slash : "";
  test_srcdir_prefix = (argc >= 3)? argv[2] + slash : test_objdir_prefix;

  try {
    test_headers(t);
    test_alignments(t);
    test_intervals(t);
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
