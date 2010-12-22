/*  samcount.cpp -- Display statistics for SAM and BAM files.

    Copyright (C) 2010 Genome Research Ltd.

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
#include <map>
#include <string>
#include <cstdlib>

#include "sam/alignment.h"
#include "sam/header.h"
#include "sam/stream.h"
#include "utilities/utilities.h"

using std::string;
using namespace sam;

struct count_pair {
  unsigned long nmapped, nunmapped;
};

typedef std::map<string, count_pair> count_map;

count_map lib_count;

string table_header(const string& div) {
  return "mapped\tunmapped\t" + div;
}

void count(isamstream& in, bool display, const string& fname) {
  in.exceptions(std::ios::failbit | std::ios::badbit);

  collection headers;
  in >> headers;

  // TODO One day we'll implement sam::read_group_header instead...
  std::map<string, string> rg_lib;
  string empty_str;
  for (collection::iterator it = headers.begin(); it != headers.end(); ++it)
    if (it->type_equals("RG"))
      rg_lib[it->field<string>("ID")] = it->field("LB", empty_str);

  alignment aln;
  string rg_buffer;

  count_map rg;

  while (in >> aln) {
    count_pair& counts = rg[aln.aux(rg_buffer, "RG", "(ungrouped)")];
    if (aln.flags() & UNMAPPED)  counts.nunmapped++;
    else  counts.nmapped++;
  }

  if (display) {
    if (! fname.empty())  std::cout << "Read groups for " << fname << ":\n";
    std::cout << table_header("readgroup") << '\n';
    for (count_map::iterator it = rg.begin(); it != rg.end(); ++it)
      std::cout << it->second.nmapped << '\t' << it->second.nunmapped << '\t'
		<< it->first << '\n';
    std::cout << "\n";
  }

  for (std::map<string, string>::iterator it = rg_lib.begin();
       it != rg_lib.end(); ++it) {
    count_pair& rg_counts = rg[it->first];
    count_pair& lib_counts = lib_count[it->second];
    lib_counts.nmapped += rg_counts.nmapped;
    lib_counts.nunmapped += rg_counts.nunmapped;
  }
}

int main(int argc, char** argv) {
  const char usage[] =
"Usage: samcount [-lr] [FILE]...\n"
"Options:\n"
"  -l  Display statistics for each library\n"
"  -r  Display statistics for each read group (by default, displays both)\n"
"";

  if (argc == 2) {
    string arg = argv[1];
    if (arg == "--version") {
      print_version(std::cout, "samcount");
      return EXIT_SUCCESS;
    }
    else if (arg == "--help") {
      std::cout << usage;
      return EXIT_SUCCESS;
    }
  }

  bool by_library = false;
  bool by_read_group = false;
  bool display_options_seen = false;

  int c;
  while ((c = getopt(argc, argv, ":lr")) >= 0)
    switch (c) {
    case 'l':  display_options_seen = by_library = true;  break;
    case 'r':  display_options_seen = by_read_group = true;  break;
    default:
      std::cerr << usage;
      return EXIT_FAILURE;
    }

  if (! display_options_seen)
    by_library = by_read_group = true;

  if (optind == argc) {
    isamstream in("-");
    count(in, by_read_group, "");
  }
  else
    for (int i = optind; i < argc; i++) {
      isamstream in(argv[i]);
      if (in.is_open())
	count(in, by_read_group, argv[i]);
      else
	std::cerr << "error opening " << argv[i] << " or something\n";
    }

  if (by_library) {
    std::cout << table_header("library") << '\n';
    for (count_map::iterator it = lib_count.begin();
	 it != lib_count.end(); ++it)
      std::cout << it->second.nmapped << '\t' << it->second.nunmapped << '\t'
		<< it->first << '\n';
  }

  return EXIT_SUCCESS;
}
