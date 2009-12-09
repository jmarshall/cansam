/** @file  sam/exception.h
    @brief Exceptions thrown by Cansam classes and functions

Blah blah.  */

#ifndef CANSAM_EXCEPTION_H
#define CANSAM_EXCEPTION_H

#include <string>
#include <ios>

namespace sam {

/** @class sam::exception sam/exception.h
    @brief Base class for Cansam exceptions */
class exception : public std::ios_base::failure {
public:
  explicit exception(const std::string& message) : failure(message) { }
  virtual ~exception() throw() { }
  //using failure::what;

  /// The filename associated with this problem, or empty if none or unknown
  std::string filename() const { return filename_; }

  /// Set an associated filename
  void set_filename(const std::string& filename) { filename_ = filename; }

private:
  std::string filename_;
};

/** @class sam::failure sam/exception.h
    @brief Exceptions representing formatting errors, etc */
class failure : public exception {
public:
  explicit failure(const std::string& message) : exception(message) { }
  virtual ~failure() throw() { }
  //using exception::what;
};

/** @class sam::sysbork sam/exception.h
    @brief Exceptions raised from system call failures */
class sysbork : public exception {
public:
  explicit sysbork(const std::string& message, int errnum)
    : exception(message), errnum_(errnum) { }
  virtual ~sysbork() throw() { }
  virtual const char* what() const throw();

  /// The @c errno error code underlying this exception
  int errnum() const { return errnum_; }

private:
  int errnum_;
};

} // namespace sam

#endif
