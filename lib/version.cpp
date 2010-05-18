/*  version.cpp
    Cansam library version information
*/

#include "sam/version.h"

#include <sstream>

namespace sam {

std::string version() {
  int major =  CANSAM_VERSION / 100000;
  int minor = (CANSAM_VERSION / 100) % 1000;
  int patch =  CANSAM_VERSION % 100;

  std::ostringstream s;
  s << major << '.' << minor;
  if (patch != 0)  s << '.' << patch;
  return s.str();
}

} // namespace sam
