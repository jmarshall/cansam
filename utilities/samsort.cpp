#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cstdlib>
#include <cstring>

#include "samsort.h"
#include "sam/alignment.h"

using sam::scoord_t;

#include "samsort.h"

using std::string;

bool lt_qname(const sam::alignment& a, const sam::alignment& b) {
  int dqname = strcmp(a.qname_c_str(), b.qname_c_str());
  if (dqname != 0)  return (dqname < 0);

  int dorder = a.order() - b.order();
  return (dorder < 0);
}

int cmp_qname(const sam::alignment& a, const sam::alignment& b) {
  int dqname = strcmp(a.qname_c_str(), b.qname_c_str());
  if (dqname != 0)  return dqname;

  int dorder = a.order() - b.order();
  return dorder;
}

bool lt_location(const sam::alignment& a, const sam::alignment& b) {
  // FIXME unmapped
  int drindex = a.rindex() - b.rindex();
  if (drindex != 0)  return (drindex < 0);

  scoord_t dpos = a.pos() - b.pos();
  if (dpos != 0)  return (dpos < 0);

  return 0;
}

// Compare by chromosome, position, read name.
int cmp_location(const sam::alignment& a, const sam::alignment& b) {
  // Treating these as unsigned means -1 (unmapped) sorts last.
  unsigned a_rindex = unsigned(a.rindex());
  unsigned b_rindex = unsigned(b.rindex());
  if (a_rindex != b_rindex)  return (a_rindex < b_rindex)? -1 : +1;

  if (a.pos() != b.pos())  return (a.pos() < b.pos())? -1 : +1;

  return cmp_qname(a, b);
}

static const alignment_comparator 
  rname_pos("location","Order by chromosome then position (and then read name)",
	    lt_location),
  qname("qname", "Order by read (query) name then first/second ordering flags",
	lt_qname);

typedef std::map<std::string, const alignment_comparator*> comparator_map;

static comparator_map& comparators() {
  static comparator_map* cmap = new comparator_map();
  return *cmap;
}

alignment_comparator::alignment_comparator(const char* name,
    const char* description0, compare* cmp0)
  : description(description0), comparer(cmp0) {
  comparators()[name] = this;
}

int main(int argc, char** argv) {
  const char usage[] =
"Usage: samsort [-bcm] [-f CMP] [-o FILE] [-S SIZE] [-T DIR] [-z NUM] [FILE]...\n"
"Options:\n"
"  -b         Write output in BAM format\n"
"  -c         Check whether input is already sorted\n"
"  -f CMP     Compare records according to comparison function CMP [location]\n"
"  -m         Merge already-sorted files\n"
"  -o FILE    Write output to FILE rather than standard output\n"
"  -S SIZE    Use SIZE amount of in-memory working space\n"
"  -T DIR     Write temporary files to DIR [$TMPDIR or /tmp]\n"
"  -z NUMBER  Compress output at level NUMBER [SAM: no compression; BAM: 6]\n"
"";

  if (argc == 2) {
    string arg = argv[1];
    if (arg == "--version") {
      std::cout << "samsort 0.1\n";
      return EXIT_SUCCESS;
    }
    else if (arg == "--help") {
      std::cout << usage;
      std::cout << "Comparison functions:\n";
      for (std::map<std::string, const alignment_comparator*>::
	       const_iterator it = comparators().begin();
	   it != comparators().end(); ++it)
	std::cout << "  " << std::left << std::setw(9) << it->first
		  << "  " << it->second->description << "\n";

      return EXIT_SUCCESS;
    }
  }

  int c;
  while ((c = getopt(argc, argv, ":bcf:mo:S:T:z:")) >= 0)
    switch (c) {
    }

  return EXIT_SUCCESS;
}
