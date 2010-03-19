#include <iostream>
#include <map>
#include <string>
#include <cstdlib>

#include "sam/alignment.h"
#include "sam/header.h"
#include "sam/stream.h"

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
      std::cout << "samcount 0.2\n";
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
