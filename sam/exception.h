/// @file sam/exception.h
/// Exceptions thrown by Cansam classes and functions

/*  Copyright (C) 2010 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
 3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
    nor the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND ITS CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */

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
  system_error(const std::string& message, int errnum)
    : exception(message), errnum_(errnum) { }

  /// Construct an exception with the given @a message, @a filename,
  /// and @a errno value
  system_error(const std::string& message, const std::string& filename,
	       int errnum)
    : exception(message), errnum_(errnum) { set_filename(filename); }

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
