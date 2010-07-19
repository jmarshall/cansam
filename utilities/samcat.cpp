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
    out << aln;
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

int main(int argc, char** argv) {
  const char usage[] =
"Usage: samcat [-bnv] [-o FILE] [FILE]...\n"
"Options:\n"
"  -b       Write output in BAM format\n"
"  -n       Suppress '@' headers in the output\n"
"  -o FILE  Write to FILE rather than standard output\n"
"  -v       Display file information and statistics\n"
"";

  string output_fname = "-";
  std::ios::openmode output_mode = sam_format;
  bool suppress_headers = false;
  bool verbose = false;

  if (argc == 1) {
    std::cerr << usage;
    return EXIT_FAILURE;
    }
  else if (argc == 2) {
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

  int c;
  while ((c = getopt(argc, argv, ":bno:v")) >= 0)
    switch (c) {
    case 'b':  output_mode = bam_format;  break;
    case 'n':  suppress_headers = true;  break;
    case 'o':  output_fname = optarg;  break;
    case 'v':  verbose = true;  break;
    default:
      std::cerr << usage;
      return EXIT_FAILURE;
    }

  try {
    osamstream out(output_fname, std::ios::out | output_mode);
    if (!out.is_open())
      throw sam::system_error("can't write to ", output_fname, errno);

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
  }
  catch (const sam::exception& e) {
    std::cerr << "samcat: " << e.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
