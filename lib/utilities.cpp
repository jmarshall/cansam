#include "lib/utilities.h"

#include <string>

using std::string;

namespace sam {

// Removes a trailing line terminator, whether it be LF, CR, or CR-LF.
// (Usually CR would be because there was a CR-LF terminator and the LF has
// already been elided.)
string& chomp(string& s) {
  string::size_type len = s.length();
  if (len >= 2 && s.compare(len - 2, 2, "\x0d\x0a") == 0)
    s.erase(len - 2);
  else if (len >= 1 && (s[len - 1] == '\x0d' || s[len - 1] == '\x0a'))
    s.erase(len - 1);

  return s;
}

coord_t to_int(const string& str,
	       string::size_type begin, string::size_type end) {
  const char* data = str.data();
  const char* s = data + begin;
  const char* lim = data + end;

  coord_t val = 0;
  while (s < lim)
    val = 10 * val + *s++ - '0';
  return val;
}

} // namespace sam
