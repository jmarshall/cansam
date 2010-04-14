#include <iostream>
#include <string>
#include <cstdlib>
#include <cerrno>

#include "sam/alignment.h"
#include "sam/exception.h"
#include "sam/header.h"
#include "sam/stream.h"

using std::string;
using namespace sam;

void group_alignments(isamstream& in, osamstream& out) {
}

int main(int argc, char** argv) {
  const char usage[] =
"Usage: samgroupbyname [-bv] [-o FILE] [FILE]\n"
"Options:\n"
"  -b       Write output in BAM format\n"
"  -o FILE  Write to FILE rather than standard output\n"
"  -v       Display file information and statistics\n"
"";

  if (argc == 2) {
    string arg = argv[1];
    if (arg == "--version") {
      std::cout << "samgroupbyname 0.3\n";
      return EXIT_SUCCESS;
    }
    else if (arg == "--help") {
      std::cout << usage;
      return EXIT_SUCCESS;
    }
  }

  string output_fname = "-";
  std::ios::openmode output_mode = sam_format;
  bool verbose = false;

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

  string input_fname = (optind < argc)? argv[optind++] : "-";
  if (optind < argc) {
    // This utility accepts at most one input file.
    std::cerr << usage;
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
  }
  catch (const sam::exception& e) {
    std::cerr << "samgroupbyname: " << e.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
