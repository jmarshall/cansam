#include <iostream>
#include <string>
#include <cstdlib>

#include "sam/alignment.h"
#include "sam/header.h"
#include "sam/stream.h"

using std::string;
using namespace sam;

void cat(isamstream& in, osamstream& out) {
  in.exceptions(std::ios::failbit | std::ios::badbit);

  collection headers;
  in >> headers;

  //std::cout << std::showpoint << headers;
  out << headers;

  alignment aln;
//  std::cout << "default ctored: "; aln.dump_on(std::cout);
  while (in >> aln) {
    //std::cout << "in loop: "; aln.dump_on(std::cout);
    out << aln;
    //std::cout << aln << '\n';
  }
}

int main(int argc, char** argv) {
  const char usage[] =
"Usage: samcat [-bv] [-o FILE] [FILE]...\n"
"Options:\n"
"  -b       Write output in BAM format\n"
"  -o FILE  Write to FILE rather than standard output\n"
"  -v       Display file information and statistics\n"
"";

  string output_fname = "-";
  std::ios::openmode output_mode = sam_format;
  bool verbose = false;

  if (argc == 1) {
    std::cerr << usage;
    return EXIT_FAILURE;
    }
  else if (argc == 2) {
    string arg = argv[1];
    if (arg == "--version") {
      std::cout << "samcat 0.1\n";
      return EXIT_SUCCESS;
    }
    else if (arg == "--help") {
      std::cout << usage;
      return EXIT_SUCCESS;
    }
  }

  int c;
  while ((c = getopt(argc, argv, ":bo:v")) >= 0)
    switch (c) {
    case 'b':  output_mode = bam_format;  break;
    case 'o':  output_fname = optarg;  break;
    case 'v':  verbose = true;  break;
    default:
      std::cerr << usage;
      return EXIT_FAILURE;
    }

  osamstream out(output_fname, std::ios::out | output_mode);

  if (optind == argc) {
    isamstream in("-");
    cat(in, out);
  }
  else
    for (int i = optind; i < argc; i++) {
      isamstream in(argv[i]);
      if (in.is_open())
	cat(in, out);
      else
	std::cerr << "error opening " << argv[i] << " or something\n";
    }

  return EXIT_SUCCESS;
}
