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
  /// Construct an exception with the given @a message
  explicit exception(const std::string& message) : failure(message) { }
  virtual ~exception() throw() { }

  //using failure::what;

  /// The filename associated with this problem, or empty if none or unknown
  std::string filename() const { return filename_; }

  /// Set an associated filename
  void set_filename(const std::string& filename) { filename_ = filename; }

protected:
  /// A buffer for use by what()
  std::string message;

private:
  std::string filename_;
};

/** @class sam::failure sam/exception.h
    @brief Exception representing formatting errors, etc */
class failure : public exception {
public:
  /// Construct a failure exception with the given @a message
  explicit failure(const std::string& message) : exception(message) { }
  virtual ~failure() throw() { }
  //using exception::what;
};

/** @class sam::system_error sam/exception.h
    @brief Exception raised from system call failures */
class system_error : public exception {
public:
  /// Construct an exception with the given @a message and @a errno value
  explicit system_error(const std::string& message, int errnum)
    : exception(message), errnum_(errnum) { }

  virtual ~system_error() throw() { }

  /// Returns a message, including filename and @a errno, if available
  virtual const char* what() const throw();

  /// The @c errno error code underlying this exception
  int errnum() const { return errnum_; }

private:
  int errnum_;
};

/** @class sam::bad_format sam/exception.h
    @brief Exception representing a SAM or BAM parsing error

Invalid records encountered while reading SAM or BAM streams cause the
stream's @c failbit state flag to be set and, if @c failbit is listed in its
@c exceptions(), a sam::bad_format exception will be thrown accordingly.  */
class bad_format : public exception {
public:
  /// Construct an exception with the given particulars
  explicit bad_format(const std::string& message, int recnum = 0)
    : exception(message), recnum_(recnum) { }

  virtual ~bad_format() throw() { }

  /// The full message associated with this problem
  virtual const char* what() const throw();

  /// The record number associated with this problem, or 0 if none or unknown
  int recnum() const { return recnum_; }

  /// Set an associated record number
  void set_recnum(int recnum) { recnum_ = recnum; }

private:
  int recnum_;
};

} // namespace sam

#endif
