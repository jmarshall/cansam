/*  samhead.cpp -- Print headers of SAM and BAM files.

    Copyright (C) 2012 Genome Research Ltd.

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
#include <string>
#include <cerrno>
#include <cstdlib>

#include "cansam/exception.h"
#include "cansam/sam/header.h"
#include "cansam/sam/stream.h"
#include "tools/utilities.h"

using std::string;
using namespace sam;

void head(const string& filename) {
  isamstream in(filename);
  if (! in.is_open())  throw sam::system_error("can't open ", filename, errno);
  in.exceptions(std::ios::failbit | std::ios::badbit);

  collection headers;
  in >> headers;

  std::cout << headers;
}

int main(int argc, char** argv)
try {
  const char usage[] =
"Usage: samhead [FILE]...\n"
"";

  if (argc >= 2) {
    string arg = argv[1];
    if (arg == "--help") { std::cout << usage; return EXIT_SUCCESS; }
    else if (arg == "--version")
      { print_version(std::cout, "samhead"); return EXIT_SUCCESS; }
  }

  if (argc == 1 && cin_likely_from_user())
    { std::cerr << usage; return EXIT_FAILURE; }

  if (argc == 1)  head("-");
  else if (argc == 2)  head(argv[1]);
  else
    for (int i = 1; i < argc; i++) {
      if (i > 1)  std::cout << '\n';
      std::cout << "==> " << argv[i] << " <==\n";
      head(argv[i]);
    }

  return EXIT_SUCCESS;
}
catch (const std::exception& e) {
  std::cout << std::flush;
  std::cerr << "samhead: " << e.what() << std::endl;
  return EXIT_FAILURE;
}
