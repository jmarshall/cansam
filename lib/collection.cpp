#include "sam/header.h"
#include "sam/exception.h"
#include "lib/sambamio.h"  // for push_back() flags

#include <iostream> // FIXME NUKE-ME

#include "lib/utilities.h"

using std::string;

namespace sam {

collection::collection() {
  allocate_cindex();
  refseqs_in_headers = false;
}

collection::~collection() {
  // delete through the pointers

  free_cindex();
}

std::vector<collection*>
  collection::collections(1, static_cast<collection*>(NULL));

void collection::allocate_cindex() {
  // TODO Reuse old ones after a while
  cindex = collections.size();
  collections.push_back(this);
}

void collection::clear() {
  headers.clear();
  refseqs.clear();
  refnames.clear();
  refseqs_in_headers = false;
}

// FIXME Make me class static?
static refsequence unmapped_refseq("*", 0, -1);

refsequence& collection::findseq(int index) {
  if (index >= 0 && index < int(refseqs.size()))
    return *refseqs[index];
  else if (index == -1)
    return unmapped_refseq;
  else
    throw sam::exception("Reference sequence index out of range"); // FIXME
}

refsequence& collection::findseq(const string& name) {
  refname_map::iterator it = refnames.find(name);
  if (it != refnames.end())
    return *(it->second);
  else if (name == "*")  // FIXME Is "*" going to be in refnames, or not?
    return unmapped_refseq;
  else
    throw sam::exception(make_string()
	<< "No such reference sequence ('" << name << "')");
}

refsequence& collection::findseq(const char* name) {
  if (name[0] == '*' && name[1] == '\0')
    return unmapped_refseq;

  return *refnames[string(name)];
}

void collection::push_back(const std::string& header_line) {
std::clog << "collection::push_back(\"" << header_line << "\")\n";
  // FIXME
}

// TEXT is NUL-delimited.
void collection::push_back(const std::string& text, int flags) {
  header* hdr;
  if (text.compare(0, 3, "@SQ") == 0) {
    refsequence* rhdr =
	new refsequence(text, (flags & add_refseq)? refseqs.size() : -1);

    if (flags & add_refseq)   refseqs.push_back(rhdr);
    if (flags & add_refname)  refnames[rhdr->name()] = rhdr;
    hdr = rhdr;
  }
  else
    hdr = new header(text);

  if (flags & add_header)  headers.push_back(hdr);
}

} // namespace sam
