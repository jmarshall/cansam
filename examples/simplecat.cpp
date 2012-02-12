#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cerrno>

#include "cansam/sam/header.h"
#include "cansam/sam/alignment.h"

void cat(std::istream& in, bool print_headers) {
  // FIXME there is no hope, as istream >> aln is silly to implement
#if 0
  sam::header header;
  while (in >> header)
    if (print_headers)
      std::cout << header << '\n';

  in.clear();

  sam::alignment aln;
  while (in >> aln)
    std::cout << aln << '\n';
#endif
}

int main(int argc, char** argv) {
  int status = EXIT_SUCCESS;

  if (argc == 1)
    cat(std::cin, true);
  else
    for (int i = 1; i < argc; i++) {
      std::ifstream strm(argv[i]);
      if (! strm.fail())
	cat(strm, (i == 1));
      else {
	std::cerr << "simplecat: can't open '" << argv[i] << "': "
		  << strerror(errno) << "\n";
	status = EXIT_FAILURE;
      }
    }

  return status;
}
