/*  samcat.cpp -- Concatenate and print SAM and BAM files.

    Copyright (C) 2010-2011 Genome Research Ltd.

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
#include <cstdlib>
#include <cerrno>

#include "sam/alignment.h"
#include "sam/exception.h"
#include "sam/header.h"
#include "sam/stream.h"
#include "utilities/utilities.h"

using std::string;
using namespace sam;

static struct samcat_options {
  int pos_flags, neg_flags;
} opt;

inline bool should_emit(const alignment& aln) {
  int flags = aln.flags();
  return (flags & opt.pos_flags) == opt.pos_flags &&
	 (flags & opt.neg_flags) == 0;
}

static struct samcat_statistics {
  unsigned long nin, nout;
} stats;

void cat(isamstream& in, osamstream& out, bool suppress_headers) {
  in.exceptions(std::ios::failbit | std::ios::badbit);

  collection headers;
  in >> headers;

  if (! suppress_headers) {
    // FIXME In BAM-world, this will have to be done with headers.clear() or so
    //std::cout << std::showpoint << headers;
    out << headers;
  }

  alignment aln;
  while (in >> aln) {
    stats.nin++;
    if (should_emit(aln)) { out << aln; stats.nout++; }
  }
}

void cat_to_fastq(isamstream& in, std::ostream& out) {
  in.exceptions(std::ios::failbit | std::ios::badbit);

  collection headers;
  in >> headers;

  alignment aln;
  string seq_buffer, qual_buffer;
  while (in >> aln) {
    out << '@' << aln.qname_c_str();
    if (aln.flags() & FIRST_IN_PAIR)  out << "/1";
    else if (aln.flags() & SECOND_IN_PAIR)  out << "/2";
    out << '\n';

    // TODO Revcomp if mapped on the negative strand

    out << aln.seq(seq_buffer) << "\n+\n" << aln.qual(qual_buffer) << '\n';
  }
}

void parse_format(const string& s,
		  std::ios::openmode& mode, std::ios::fmtflags& format) {
  if (s == "bam")  mode = bam_format;
  else if (s == "hex")  format = std::ios::hex;
  else if (s == "text")  format = std::ios::boolalpha;
  else  throw bad_format("Invalid output format ('" + s + "')");
}

int main(int argc, char** argv)
try {
  const char usage[] =
"Usage: samcat [-bnv] [-f FLAGS] [-o FILE] [-O FORMAT] [FILE]...\n"
"Options:\n"
"  -b         Write output in BAM format (equivalent to -Obam)\n"
"  -f FLAGS   Display only alignment records matching FLAGS\n"
"  -n         Suppress '@' headers in the output\n"
"  -o FILE    Write to FILE rather than standard output\n"
"  -O FORMAT  Write output in the specified FORMAT\n"
"  -v         Display file information and statistics\n"
"Output formats:\n"
"  bam        Compressed binary BAM format\n"
"  hex        SAM format, with flags displayed in hexadecimal\n"
"  text       SAM format, with flags displayed as readable strings\n"
"";

  string output_fname = "-";
  std::ios::openmode output_mode = sam_format;
  std::ios::fmtflags output_format = std::ios::dec;
  bool suppress_headers = false;
  bool verbose = false;

  if (argc == 2) {
    string arg = argv[1];
    if (arg == "--version") {
      print_version(std::cout, "samcat");
      return EXIT_SUCCESS;
    }
    else if (arg == "--help") {
      std::cout << usage;
      return EXIT_SUCCESS;
    }
  }

  opt.pos_flags = opt.neg_flags = 0;

  int c;
  while ((c = getopt(argc, argv, ":bf:no:O:v")) >= 0)
    switch (c) {
    case 'b':  output_mode = bam_format;  break;
    case 'f':  parse_flags(optarg, opt.pos_flags, opt.neg_flags);  break;
    case 'n':  suppress_headers = true;  break;
    case 'o':  output_fname = optarg;  break;
    case 'O':  parse_format(optarg, output_mode, output_format);  break;
    case 'v':  verbose = true;  break;
    default:
      std::cerr << usage;
      return EXIT_FAILURE;
    }

  if (argc == 1 && cin_likely_from_user())
    { std::cerr << usage; return EXIT_FAILURE; }

  stats.nin = stats.nout = 0;

  osamstream out(output_fname, std::ios::out | output_mode);
  if (!out.is_open())
    throw sam::system_error("can't write to ", output_fname, errno);

  out.setf(output_format, std::ios::basefield | std::ios::boolalpha);

  if (optind == argc) {
    isamstream in("-");
    cat(in, out, suppress_headers);
  }
  else
    for (int i = optind; i < argc; i++) {
      isamstream in(argv[i]);
      if (in.is_open())
	cat(in, out, suppress_headers);
      else {
	sam::system_error error("can't open ", argv[i], errno);
	std::cerr << "samcat: " << error.what() << '\n';
      }
    }

  if (verbose) {
    std::clog << "Wrote " << stats.nout << " records";
    if (stats.nout != stats.nin)  std::clog << " (out of " << stats.nin << ')';
    if (output_fname != "-")  std::clog << " to " << output_fname;
    std::clog << '\n';
  }

  return EXIT_SUCCESS;
}
catch (const sam::exception& e) {
  std::cout << std::flush;
  std::cerr << "samcat: " << e.what() << std::endl;
  return EXIT_FAILURE;
}
