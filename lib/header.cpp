#include "sam/header.h"

#include <string>

using std::string;

namespace sam {

header& header::assign(const string& line) {
  // FIXME hmmm?
  // Accept a missing '@' to avoid inventing an error to be handled.
  string::size_type pos = (line.length() >= 1 && line[0] == '@')? 1 : 0;

  string::size_type tabpos = line.find('\t', pos);
  type_.assign(line, pos, tabpos - pos);

  clear();
  while (tabpos != string::npos) {
    pos = tabpos + 1;
    tabpos = line.find('\t', pos);

    if (tabpos - pos >= 3 && line[pos + 2] == ':')
      push_back(header_field(line.substr(pos, 2),
			     line.substr(pos + 3, tabpos - (pos + 3))));
    else {
      // The field has no tag; sadly this can happen, e.g., in @CO comments.
      push_back("", line.substr(pos, tabpos - pos));
    }
  }

  return *this;
}

string header::str() const {
  string s = "@";
  s += type_;

  for (const_iterator it = begin(); it != end(); ++it) {
    s += '\t';
    if (! it->tag.empty()) {
      s += it->tag;
      s += ':';
      }
    s += it->value;
  }

  return s;
}

header::const_iterator header::find(const char* tag) const {
  for (const_iterator it = begin(); it != end(); ++it)
    if (it->tag == tag)  return it;

  return end();
}

header::iterator header::find(const char* tag) {
  for (iterator it = begin(); it != end(); ++it)
    if (it->tag == tag)  return it;

  return end();
}

std::string header::field(const char* tag,
			  const std::string& default_value) const {
  for (const_iterator it = begin(); it != end(); ++it)
    if (it->tag == tag)  return it->value;

  return default_value;
}

} // namespace sam
