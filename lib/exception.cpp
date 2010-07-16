#include "sam/exception.h"

#include <sstream>
#include <cstring>  // for strerror()

namespace sam {

namespace {

// Returns whether S ends in a non-space character, i.e., is a complete phrase.
bool complete_phrase(const std::string& s) {
  return (! s.empty()) && s[s.length() - 1] != ' ';
}

} // unnamed namespace

// Returns "message [[for ] "filename"] [[at record] N]"
const char* bad_format::what() const throw() {
  std::ostringstream s;

  s << sam::exception::what();
  if (! filename().empty()) {
    if (complete_phrase(s.str()))  s << " for ";
    s << '"' << filename() << '"';
  }

  if (recnum() != 0) {
    if (complete_phrase(s.str()))  s << " at record ";
    s << recnum();
  }

  what_text_ = s.str();
  return what_text_.c_str();
}

// Returns "message [[for ] "filename"]: strerror
const char* system_error::what() const throw() {
  what_text_ = sam::exception::what();

  if (! filename().empty()) {
    if (complete_phrase(what_text_))  what_text_ += " for ";
    what_text_ += '"';
    what_text_ += filename();
    what_text_ += '"';
  }

  what_text_ += ": ";
  what_text_ += strerror(errnum());

  return what_text_.c_str();
}

} // namespace sam
