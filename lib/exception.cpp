#include "sam/exception.h"

namespace sam {

const char* sysbork::what() const throw() {
  return "foobar";
}

} // namespace sam
