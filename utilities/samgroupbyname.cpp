/*  samgroupbyname.cpp -- Order a SAM/BAM file so that read pairs are together.

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

#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <cerrno>
#include <cstdlib>

#include "sam/algorithm.h"
#include "sam/alignment.h"
#include "sam/exception.h"
#include "sam/header.h"
#include "sam/stream.h"
#include "utilities/utilities.h"

using std::string;
using namespace sam;

typedef std::set<alignment, less_by_qname> alignment_set;

struct group_statistics {
  group_statistics() { pairs = singletons = max_pending = 0; }

  unsigned long pairs;
  unsigned long singletons;
  unsigned long max_pending;
};

bool emit_singletons = true;
bool verbose = false;
group_statistics stats;

void group_alignments(isamstream& in, osamstream& out) {
  alignment_set seen;
  unsigned long seen_size = 0;

  alignment aln;
  while (in >> aln) {
    std::pair<alignment_set::iterator, bool> ins = seen.insert(aln);
    if (ins.second) {
      seen_size++;
    }
    else {
      stats.pairs++;
      out << *ins.first << aln;
      seen.erase(ins.first);

      if (seen_size > stats.max_pending)  stats.max_pending = seen_size;
      seen_size--;
    }
  }

  stats.singletons += seen_size;
  if (emit_singletons)
    for (alignment_set::iterator it = seen.begin(); it != seen.end(); ++it)
      out << *it;
}

int main(int argc, char** argv) {
  const char usage[] =
"Usage: samgroupbyname [-bpv] [-o FILE] [FILE]\n"
"Options:\n"
"  -b       Write output in BAM format\n"
"  -o FILE  Write to FILE rather than standard output\n"
"  -p       Emit pairs only, discarding any leftover singleton reads\n"
"  -v       Display file information and statistics\n"
"";

  if (argc == 2) {
    string arg = argv[1];
    if (arg == "--version") {
      print_version(std::cout, "samgroupbyname");
      return EXIT_SUCCESS;
    }
    else if (arg == "--help") {
      std::cout << usage;
      return EXIT_SUCCESS;
    }
  }

  string output_fname = "-";
  std::ios::openmode output_mode = sam_format;

  int c;
  while ((c = getopt(argc, argv, ":bo:pv")) >= 0)
    switch (c) {
    case 'b':  output_mode = bam_format;  break;
    case 'o':  output_fname = optarg;  break;
    case 'p':  emit_singletons = false;  break;
    case 'v':  verbose = true;  break;
    default:
      std::cerr << usage;
      return EXIT_FAILURE;
    }

  string input_fname = (optind < argc)? argv[optind++] : "-";
  if (optind < argc) {
    std::cerr << "samgroupbyname: "
		 "only one input file can be processed at a time" << std::endl;
    return EXIT_FAILURE;
  }

  try {
    isamstream in(input_fname);
    if (!in.is_open())
      throw sam::system_error("can't open ", input_fname, errno);

    osamstream out(output_fname, output_mode);
    if (!out.is_open())
      throw sam::system_error("can't write to ", output_fname, errno);

    in.exceptions (std::ios::failbit | std::ios::badbit);
    out.exceptions(std::ios::failbit | std::ios::badbit);

    collection headers;
    in >> headers;

    // After grouping by read name, the file will no longer be sorted;
    // so replace any sort-order tags by the appropriate group-order tag.
    for (collection::iterator it = headers.begin(); it != headers.end(); ++it)
      if (it->type_equals("HD")) {
	it->erase("SO");
	it->set_field("GO", "query");
      }

    out << headers;
    group_alignments(in, out);

    if (verbose) {
      const char* action = emit_singletons? "written:   " : "discarded: ";
      std::clog
	<<   "Paired reads written:     " << std::setw(12) << stats.pairs * 2
	<< "\nUnpaired reads " << action  << std::setw(12) << stats.singletons
	<< "\nMaximum reads in memory:  " << std::setw(12) << stats.max_pending
	<< std::endl;
    }
  }
  catch (const sam::exception& e) {
    std::cerr << "samgroupbyname: " << e.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
