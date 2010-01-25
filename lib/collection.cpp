#include "sam/header.h"
#include "sam/collection.h"
#include "sam/exception.h"

using std::string;

namespace sam {

collection::collection() {
}

collection::~collection() {
  // delete through the pointers

  collections[cindex] = NULL;
}

std::vector<collection*> collection::collections;

void collection::clear() {
  // FIXME
}

// FIXME Make me class static?
static refsequence unmapped_refseq("*", 0, -1);

refsequence& collection::findseq(int index) {
  if (index >= 0 && size_t(index) < refseqs.size())
    return *refseqs[index];
  else if (index == -1)
    return unmapped_refseq;
  else
    throw bad_format("Reference sequence index out of range"); // FIXME
}

void collection::push_back(const std::string& header_line) {
  // FIXME
}

} // namespace sam
