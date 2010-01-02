#include "sam/exception.h"

namespace sam {

const char* system_error::what() const throw() {
  return "FIXME system_error";
}

const char* bad_format::what() const throw() {
  return std::ios_base::failure::what();
  return "FIXME bad_format";
}

} // namespace sam
