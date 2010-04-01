#include "sam/exception.h"

#include <cstring>  // for strerror()

namespace sam {

const char* bad_format::what() const throw() {
  // FIXME
  what_text_ = sam::exception::what();
  what_text_ += " (bad_format)";
  return what_text_.c_str();
}

// Returns "message [[for ] "filename"]: strerror
const char* system_error::what() const throw() {
  what_text_ = sam::exception::what();

  if (! filename().empty()) {
    if (! (what_text_.empty() || what_text_[what_text_.length() - 1] == ' '))
      what_text_ += " for ";
    what_text_ += '"';
    what_text_ += filename();
    what_text_ += '"';
  }

  what_text_ += ": ";
  what_text_ += strerror(errnum());

  return what_text_.c_str();
}

} // namespace sam
