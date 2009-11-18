#include <iostream>
#include <string>
#include <cstdlib>

#include "sam/header.h"
#include "sam/alignment.h"

using std::string;

void cat(std::istream& in, std::ostream& out) {
  sam::header header;
  while (in >> header)
    out << header << '\n';

  in.clear();

  sam::alignment aln;
  while (in >> aln) {
    out << aln << '\n';
#if 0
    std::clog << "name{" << aln.qname << "} flags{" << aln.flags << "} rname{"
	<< aln.rname << "} pos{" << aln.pos << "} mapq{" << aln.mapq
	<< "} cigar{" << aln.cigar << "} mrname{" << aln.mate_rname
	<< "} mpos{" << aln.mate_pos << "} isize{" << aln.isize
	<< "} seq{" << aln.seq << "} qual{" << aln.qual << "}";
    if (! aln.extras.empty())
	std::clog << " extras{" << aln.extras << "}";
    std::clog << '\n';
#endif
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

  bool bam_output = false;
  string output_fname;
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
    case 'b':  bam_output = true;  break;
    case 'o':  output_fname = optarg;  break;
    case 'v':  verbose = true;  break;
    default:
      std::cerr << usage;
      return EXIT_FAILURE;
    }

  if (optind == argc)
    cat(std::cin, std::cout);
  else
    for (int i = optind; i < argc; i++) {
      string arg = argv[i];
      if (arg == "-")
	cat(std::cin, std::cout);
      else {
	//std::ifstream str
	std::cout << "processing {" << arg << "}\n";
      }
    }

  return EXIT_SUCCESS;
}
