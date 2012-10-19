/*  test/sam.cpp -- Tests for SAM and BAM formatting.

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

#include <iostream>
#include <sstream>

#include "cansam/sam/alignment.h"
#include "cansam/sam/stream.h"
#include "test/test.h"

static void test_reader(test_harness& t) {
  std::istringstream sam1(
"foo\t37\t*\t0\t0\t*\t*\t0\t0\tATGC\t????\tNM:i:4\n");

  sam::isamstream str(sam1.rdbuf());
  str.exceptions(std::ios::failbit | std::ios::badbit);
  sam::alignment aln;
std::cout << "loop 1\n";
  while (str >> aln) {
    aln.dump_on(std::cout);
    std::cout << aln << '\n';
  }
std::cout << "end of loop 1\n";

  std::cout << "* from /dev/null:\n";
  sam::isamstream str2("/dev/null");
  while (str2 >> aln)
    std::cout << aln << '\n';
}

static void test_bam_headers(test_harness& t, const string& basename,
			     const std::stringstream& text) {
  string filename = test_objdir_prefix + basename + "-out.bam";

  sam::isamstream sam(text.rdbuf());
  sam::collection headers;
  sam >> headers;

  sam::osamstream out(filename, sam::bam_format);
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out << headers;
  out.close();

  sam::isamstream in(filename);
  in.exceptions(std::ios::failbit | std::ios::badbit);

  sam::collection headers2;
  in >> headers2;
  std::ostringstream text2;
  text2 << headers2;

  t.check(text.str(), text2.str(), basename);
}

static void check_headers(test_harness& t, sam::isamstream& in,
			  const std::vector<string>& hdtext) {
  t.check(in.is_open(), "opening " + in.filename());

  sam::collection headers;
  t.check(in >> headers, "reading headers from " + in.filename());

  bool ok = headers.size() == hdtext.size();

  sam::collection::iterator it = headers.begin();
  for (size_t i = 0; ok && i < hdtext.size(); ++it, ++i)
    ok = it->str() == hdtext[i];

  t.check(ok, "headers from " + in.filename());
}

static void test_streams(test_harness& t, const string& extension,
			 std::ios::openmode mode) {
  std::vector<string> hdtext;
  hdtext.push_back("@SQ\tSN:1\tLN:249250621");
  hdtext.push_back("@SQ\tSN:MT\tLN:16569");

  std::stringstream text;
  for (std::vector<string>::iterator it = hdtext.begin();
       it != hdtext.end(); ++it)
    text << *it << '\n';

  sam::isamstream str(text.rdbuf());
  sam::collection headers;
  str >> headers;
  if (! str)  throw "extracting headers failed";

  {
    sam::osamstream out;
    out.open(test_objdir_prefix + "flushing-close-out" + extension, mode);
    t.check(out << headers, "write to " + out.filename());
    out.close();

    out.open(test_objdir_prefix + "flushing-dtor-out" + extension, mode);
    out << headers;
    t.check(out, "write to " + out.filename());
  }

  sam::osamstream out;
  out.open(test_objdir_prefix + "flushing-flush-out" + extension, mode);
  t.check(out << headers, "write to " + out.filename());
  t.check(out.flush(), "flush " + out.filename());
  // Just flushed, not closed or destructed.

  sam::isamstream in(test_objdir_prefix + "flushing-flush-out" + extension);
  check_headers(t, in, hdtext);
  in.close();

  in.open(test_objdir_prefix + "flushing-dtor-out" + extension);
  check_headers(t, in, hdtext);
  in.close();

  in.open(test_objdir_prefix + "flushing-close-out" + extension);
  check_headers(t, in, hdtext);
}

void test_sam_io(test_harness& t) {
  test_reader(t);

  std::stringstream text;
  for (int i = 1; i <= 20000; i++)
    text << "@SQ\tSN:chr" << i << "\tLN:" << 10000 + i * 2 << '\n';
  test_bam_headers(t, "manysmallheaders", text);

  text.str("");
  for (int i = 1; i <= 1200; i++) {
    text << "@SQ\tSN:chr" << i << "\tLN:" << 10000 + i * 2 << '\n';
    if (i == 1000)
      text << "@SQ\tSN:chrBig\tLN:80000\tSP:" << string(80000, 'A') << '\n';
  }
  test_bam_headers(t, "hugeheader", text);

  test_streams(t, ".sam", sam::sam_format);
  test_streams(t, ".bam", sam::bam_format);
}
