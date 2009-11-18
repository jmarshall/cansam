/** @file  sam/exception.h
    @brief Exception classes
*/

#ifndef CANSAM_EXCEPTION_H
#define CANSAM_EXCEPTION_H

#include <exception>
#include <string>

namespace sam {

/** @brief Exceptions thrown by Cansam functions.  */
class failure : public std::exception {
public:
  explicit failure(const std::string& message) : what_(message) { }
  virtual ~failure() { }
  virtual const char* what() const throw() { return what_.c_str(); }

private:
  std::string what_;
};

} // namespace sam

#endif
