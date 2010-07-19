#include "sam/header.h"
#include "sam/exception.h"
#include "lib/sambamio.h"  // for push_back() flags

#include "lib/utilities.h"

using std::string;

namespace sam {

collection::collection() {
  allocate_cindex();
  refseqs_in_headers = false;
}

collection::~collection() {
  clear();
  free_cindex();
}

std::vector<collection*>
  collection::collections(1, static_cast<collection*>(NULL));

void collection::allocate_cindex() {
  // TODO Reuse old ones after a while
  cindex = collections.size();
  collections.push_back(this);
}

template <typename InputIterator>
void delete_each(InputIterator first, InputIterator last) {
  for (InputIterator it = first; it != last; ++it)
    delete *it;
}

void collection::clear() {
  refnames.clear();
  rgroups.clear();

  if (! refseqs_in_headers)
    delete_each(refseqs.begin(), refseqs.end());
  refseqs.clear();

  delete_each(headers.begin(), headers.end());
  headers.clear();

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

readgroup& collection::findgroup_(const std::string& id) const {
  readgroup_map::const_iterator it = rgroups.find(id);
  if (it == rgroups.end())
    throw sam::exception(make_string()
	<< "No such read group ('" << id << "')");

  return *(it->second);
}

void collection::push_back(const std::string& header_line) {
  std::string text = header_line;
  size_t pos = 0;
  while ((pos = text.find('\t', pos)) != string::npos)
    text[pos++] = '\0';

  // FIXME  Is all-flags-always correct?
  push_back(text, add_header | add_refseq | add_refname);
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
  else if (text.compare(0, 3, "@RG") == 0) {
    readgroup* rghdr = new readgroup(text);
    rgroups[rghdr->id()] = rghdr;
    hdr = rghdr;
  }
  else
    hdr = new header(text);

  if (flags & add_header)  headers.push_back(hdr);
}

} // namespace sam
