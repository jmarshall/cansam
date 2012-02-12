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
  sam::isamstream str2("/dev/null", sam::sam_format);
  while (str2 >> aln)
    std::cout << aln << '\n';
}

void test_sam_io(test_harness& t) {
  test_reader(t);
}
