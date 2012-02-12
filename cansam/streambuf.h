/// @file cansam/streambuf.h
/// Low-level classes for input/output

/*  Copyright (C) 2010-2012 Genome Research Ltd.

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

#ifndef CANSAM_STREAMBUF_H
#define CANSAM_STREAMBUF_H

#include <streambuf>

/** @file
Provides low-level input/output classes derived from @c std::streambuf.
Most code will not need to use these classes directly.  */

namespace sam {

/** @class sam::streambuf cansam/streambuf.h
    @brief Stream buffer that can be open or closed */
class streambuf : public std::streambuf {
public:
  virtual ~streambuf() { }

  /// Returns whether the underlying file has been successfully opened
  virtual bool is_open() const = 0;

  /// Close the underlying file (if it is open)
  virtual void close() = 0;

protected:
  streambuf() { }

private:
  // Prevent copy construction and assignment
  streambuf(const streambuf&);
  streambuf& operator= (const streambuf&);
};

/** @class sam::rawfilebuf cansam/streambuf.h
    @brief Unbuffered file descriptor stream buffer

Provides unbuffered bulk access to a Unix-style file descriptor.
The usual public @c std::streambuf methods are available as follows:
  - @c sgetn() corresponds directly to a @c read(2) system call;
  - @c sputn() corresponds to a @c write(2) system call, invoking it repeatedly
    if necessary until all data has been written;
  - @c pubseekoff() and @c pubseekpos() correspond directly to @c lseek(2)
    system calls;
  - @c in_avail() corresponds (for ordinary files) to a combination of
    @c lseek(2) and @c fstat(2) system calls;
  - methods for character-orientated input/output (@c sgetc(), @c sputc(), etc)
    are not implemented and simply throw @c std::logic_error;
  - other methods (e.g., @c pubimbue()) have no effect.

These methods retry their system calls if they are interrupted by signal
delivery.  Thus calling code does not need to deal with @c EINTR or
foreshortened interrupted writes itself.  If a system call fails for other
reasons, these methods throw sam::system_error exceptions accordingly.

(To be precise, this class overrides the @c seekoff(), @c seekpos(),
@c showmanyc(), @c xsgetn(), and @c xsputn() protected methods, and trivially
overrides @c overflow(), @c uflow(), and @c underflow() as @c throw statements;
all the others are inherited as no-ops from @c std::streambuf.)  */
class rawfilebuf : public sam::streambuf {
public:
  /// Construct a closed buffer
  rawfilebuf() : fd_(-1), owned_(false) { }

  /// Construct a buffer, opening a file that will be closed when this buffer
  /// is destroyed
  rawfilebuf(const char* fname, std::ios_base::openmode mode, int perm = 0664)
    : fd_(-1), owned_(true) { open(fname, mode, perm); }

  /// Destroy this buffer object, optionally closing the underlying
  /// file descriptor
  ~rawfilebuf() { if (owned_)  close_nothrow(); }

  /// Open a file that will be closed when this buffer is destroyed
  rawfilebuf* open(const char* fname, std::ios_base::openmode mode,
		   int perm = 0664);

  /// Open a file that will be closed when this buffer is destroyed
  /** @param flags  flags as used by the @c open(2) system call
  */
  rawfilebuf* open(const char* fname, int flags, int perm = 0664);

  /// Associate an open file descriptor that will be closed when this buffer
  /// is destroyed
  rawfilebuf* open(int fd)
    { if (is_open())  return NULL;
      fd_ = fd; owned_ = true; return this; }

  /// Attach an open file descriptor that will not be automatically closed
  rawfilebuf* attach(int fd)
    { if (is_open())  return NULL;
      fd_ = fd; owned_ = false; return this; }

  /// Returns whether the file has been successfully opened
  virtual bool is_open() const { return fd_ >= 0; }

  /// Close the underlying file descriptor (if it is open)
  virtual void close();

  /// Returns the underlying file descriptor
  int fd() const { return fd_; }

protected:
  // @cond infrastructure
  virtual std::streamsize xsgetn(char*, std::streamsize);
  virtual std::streamsize xsputn(const char*, std::streamsize);

  virtual std::streamsize showmanyc();

  virtual std::streampos seekoff(std::streamoff, std::ios_base::seekdir,
	    std::ios_base::openmode = std::ios_base::in | std::ios_base::out);
  virtual std::streampos seekpos(std::streampos,
	    std::ios_base::openmode = std::ios_base::in | std::ios_base::out);

  virtual int_type uflow();
  virtual int_type underflow();
  virtual int_type overflow(int_type c = traits_type::eof());
  // @endcond

private:
  int fd_;
  bool owned_;

  int close_nothrow();

  // Prevent copy construction and assignment
  rawfilebuf(const rawfilebuf&);
  rawfilebuf& operator= (const rawfilebuf&);
};

} // namespace sam

#endif
