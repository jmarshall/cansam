/** @file sam/iobuffer.h
   @brief Low-level classes for input/output

foobar baz.
*/

#ifndef CANSAM_IOBUFFER_H
#define CANSAM_IOBUFFER_H

#include <ios>
#include <streambuf>

namespace sam {

/** @class sam::rawfilebuf sam/iobuffer.h
    @brief Unbuffered file descriptor stream buffer

The seekoff(), seekpos(), showmanyc(), xsgetn(), and xsputn() methods are
overridden; all the others are inherited as no-ops from std::streambuf.
*/
class rawfilebuf : public std::streambuf {
public:
  // FIXME Maybe extra parameter bool autoclose?
  rawfilebuf(int fd) : fd_(fd) { }
  ~rawfilebuf() { }

  /// Returns the underlying file descriptor
  int fd() { return fd_; }

protected:
  // FIXME Virtual functions that I maybe should be overriding:
  // imbue (default does nothing, ok)
  // overflow (for output; default does nothing, ok)
  // pbackfail (default just fails, ok)
  // setbuf (default does nothing, ok)
  // sync (default does nothing, ok)
  // uflow (default calls underflow(), ok if that returns EOF)
  // underflow (default just returns EOF, ok)

  // FIXME Don't really need to override these two
  virtual int_type uflow();
  virtual int_type underflow();

  virtual std::streamsize xsgetn(char*, std::streamsize);
  virtual std::streamsize xsputn(const char*, std::streamsize);

  virtual std::streamsize showmanyc();

  virtual std::streampos seekoff(std::streamoff, std::ios_base::seekdir,
	    std::ios_base::openmode = std::ios_base::in | std::ios_base::out);
  virtual std::streampos seekpos(std::streampos,
	    std::ios_base::openmode = std::ios_base::in | std::ios_base::out);

private:
  int fd_;
};

} // namespace sam

#endif
