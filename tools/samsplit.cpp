/*  samsplit.cpp -- Split a SAM or BAM file into separate read groups.

    Copyright (C) 2011-2012 Genome Research Ltd.

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
#include <cerrno>
#include <cstdlib>

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"
#include "cansam/sam/stream.h"
#include "cansam/exception.h"
#include "lib/utilities.h"
#include "tools/utilities.h"

using std::string;
using namespace sam;

string input_basename;
string output_extension;

string expand(const string& templ, const header& rg, int rg_index) {
  string s = templ;
  size_t pos = 0;

  while ((pos = s.find('%', pos)) != string::npos) {
    if (pos+1 >= s.length())
      throw sam::bad_format("template has orphaned '%' at end");

    int key_size = 2;
    string value;

    switch (s[pos+1]) {
    case '*':  value = input_basename;  break;
    case '.':  value = output_extension;  break;
    case '%':  value = "%";  break;

    case '#': {
      char buffer[format::buffer<int>::size];
      value.assign(buffer, format::decimal(buffer, rg_index) - buffer);
      }
      break;

    default:
      if (pos+2 >= s.length())
	throw sam::bad_format(make_string()
	    << "template ends with invalid expansion ('%" << s[pos+1] << "')");

      key_size = 3;
      value = rg.field<string>(&s[pos+1]);
      break;
    }

    s.replace(pos, key_size, value);
    pos += value.length();
  }

  return s;
}

struct split_options {
  split_options() { }

  int pos_flags, neg_flags;
  int min_quality;
};

static split_options opt;

static struct split_statistics {
  unsigned long total, discarded;
} stats = { 0, 0 };

struct split {
  split() : count(0), out(NULL) { }

  unsigned long count;
  osamstream* out;
};

typedef std::map<string, split> rg_split_map;

void split_reads(isamstream& in, rg_split_map& rg_split, osamstream* out) {
  alignment aln;
  string rg_buffer;

  while (in >> aln) {
    stats.total++;
    if ((aln.flags() & opt.pos_flags) == opt.pos_flags &&
	(aln.flags() & opt.neg_flags) == 0 &&
	aln.mapq() >= opt.min_quality) {
      split& split = rg_split[aln.aux(rg_buffer, "RG")];
      *split.out << aln;
      split.count++;

      if (out)  *out << aln;
    }
    else
      stats.discarded++;
  }
}

int main(int argc, char** argv)
try {
  static const char usage[] =
"Usage: samsplit [OPTION]... FILE [TEMPLATE]\n"
"Options:\n"
"  -b        Write output files in BAM format\n"
"  -f FLAGS  Emit only alignment records matching FLAGS\n"
"  -o FILE   Write all selected records to FILE, in addition to splitting\n"
"  -q NUM    Discard reads with mapping quality less than NUM\n"
"  -z NUM    Compress output files at level NUM (default for BAM; none for SAM)\n"
"Template and output file expansions:\n"
"  %XY       Read group header's XY field\n"
"  %#        Index of the read group (within the @RG headers, from 1)\n"
"  %*        Input FILE basename, without directory part or extension\n"
"  %.        \"sam\" or \"bam\", as appropriate for the chosen output format\n"
"  %%        A single \"%\" character\n"
"The output TEMPLATE defaults to \"%*-%ID.%.\"\n"
"";

  if (argc >= 2) {
    string arg = argv[1];
    if (arg == "--help") { std::cout << usage; return EXIT_SUCCESS; }
    else if (arg == "--version")
      { print_version(std::cout, "samsplit"); return EXIT_SUCCESS; }
  }

  string split_template = "%*-%ID.%.";
  string output_filename;

  std::ios::openmode output_mode = sam_format;
  output_extension = "sam";

  int c;
  while ((c = getopt(argc, argv, ":bf:o:q:z:")) >= 0)
    switch (c) {
    case 'b':  output_mode |= std::ios::binary; output_extension = "bam"; break;
    case 'f':  parse_flags(optarg, opt.pos_flags, opt.neg_flags);  break;
    case 'o':  output_filename = optarg;  break;
    case 'q':  opt.min_quality = atoi(optarg);  break;
    case 'z':  if (atoi(optarg) > 0)  output_mode |= compressed;
	       else  output_mode &= ~compressed;
	       break;
    default:   std::cerr << usage; return EXIT_FAILURE;
    }

  int nargs = argc - optind;
  if ((argc == 1 && cin_likely_from_user()) || nargs > 2)
    { std::cerr << usage; return EXIT_FAILURE; }

  string filename = (nargs >= 1)? argv[optind] : "-";
  if (nargs >= 2)  split_template = argv[optind+1];

  isamstream in(filename);
  if (! in.is_open())
    throw sam::system_error("can't open ", filename, errno);

  in.exceptions(std::ios::failbit | std::ios::badbit);
  input_basename = (filename != "-")? basename(filename) : "stdin";

  collection headers;
  in >> headers;

  osamstream* copy_out = NULL;
  if (! output_filename.empty()) {
    header empty("@RG");
    string copyname = expand(output_filename, empty, 0);
    copy_out = new osamstream(copyname, output_mode);
    if (! copy_out->is_open()) {
      int errno0 = errno;
      delete copy_out;
      throw sam::system_error("can't write to ", copyname, errno0);
    }

    *copy_out << headers;
  }

  // FIXME This assumes we have enough file descriptors for all the read groups.
  rg_split_map rg_split;
  int rg_index = 0;
  for (collection::iterator it = headers.begin(); it != headers.end(); ++it)
    if (it->type_equals("RG")) {
      rg_index++;
      string splitname = expand(split_template, *it, rg_index);
      osamstream* out = new osamstream(splitname, output_mode);
      if (! out->is_open()) {
	// FIXME Leaks all those allocated so far
	throw sam::system_error("can't write to ", splitname, errno);
      }
      out->exceptions(std::ios::failbit | std::ios::badbit);
      rg_split[it->field<string>("ID")].out = out;

      // TODO Remove the other @RG headers.
      *out << headers;
    }

  split_reads(in, rg_split, copy_out);

  for (rg_split_map::iterator it = rg_split.begin(); it != rg_split.end(); ++it)
    delete it->second.out;

  delete copy_out;

  return EXIT_SUCCESS;
}
catch (const std::exception& e) {
  std::cout << std::flush;
  std::cerr << "samsplit: " << e.what() << std::endl;
  return EXIT_FAILURE;
}
