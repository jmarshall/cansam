#include <string>

#include "sam/tagfield.h"

#include "lib/utilities.h"

using std::string;

namespace sam {

tagfield::tagfield(const string& s) {
  if (! (s.length() >= 3 && s[2] == ':'))
    throw "Too short or malformed tagfield"; // FIXME throw what?

  tag[0] = s[0];
  tag[1] = s[1];
  tag[2] = '\0';

  type = (s.length() >= 5 && s[4] == ':')? s[3] : '@';

  switch (type) {
  case '@':
    str.assign(s, 3, string::npos);
    break;

  case 'A':
  case 'H':
  case 'Z':
    str.assign(s, 5, string::npos);
    break;

  case 'c':
  case 'C':
  case 's':
  case 'S':
  case 'i':
  case 'I':
    i = to_int(s, 5, s.length());
    break;

  case 'f':
    abort(); // FIXME
    break;

  case 'd':
    abort(); // FIXME
    break;
  }
}

} // namespace sam
