/** @file  sam/exception.h
    @brief Exceptions thrown by Cansam classes and functions

Blah blah.  */

#ifndef CANSAM_EXCEPTION_H
#define CANSAM_EXCEPTION_H

#include <string>
#include <stdexcept>

namespace sam {

/** @class sam::exception sam/exception.h
    @brief Base class for Cansam exceptions */
class exception : public std::runtime_error {
public:
  /// Construct an exception with the given @a message
  explicit exception(const std::string& message)
    : std::runtime_error(message) { }

  virtual ~exception() throw() { }

  /// The filename associated with this problem, or empty if none or unknown
  std::string filename() const { return filename_; }

  /// Set an associated filename
  void set_filename(const std::string& filename) { filename_ = filename; }

private:
  std::string filename_;
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
  mutable std::string what_text_;
};

/** @class sam::system_error sam/exception.h
    @brief Exception raised from system call failures

The message returned by what() is in the format
@verbatim
message [[for] "filename"]: system error text
@endverbatim
where @a message is as passed to the constructor, the filename portion appears
if there is an associated filename (with "for" unless @a message ends with a
space), and the system error text is as returned by @c strerror(3).

@note Because the standard library's @c errno is a macro, this class's member
function can't be named @c errno().  */
class system_error : public exception {
public:
  /// Construct an exception with the given @a message and @a errno value
  explicit system_error(const std::string& message, int errnum)
    : exception(message), errnum_(errnum) { }

  virtual ~system_error() throw() { }

  /// Returns a message, including filename (if available) and @a errno
  virtual const char* what() const throw();

  /// The @c errno error code underlying this exception
  int errnum() const { return errnum_; }

private:
  int errnum_;
  mutable std::string what_text_;
};

} // namespace sam

#endif
