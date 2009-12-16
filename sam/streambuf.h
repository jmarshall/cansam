/** @file sam/streambuf.h
   @brief Low-level classes for input/output

foobar baz.
*/

#ifndef CANSAM_STREAMBUF_H
#define CANSAM_STREAMBUF_H

#include <streambuf>

namespace sam {

/** @class sam::streambuf sam/streambuf.h
    @brief Stream buffer that can be open or closed  */

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

/** @class sam::rawfilebuf sam/streambuf.h
    @brief Unbuffered file descriptor stream buffer

The seekoff(), seekpos(), showmanyc(), xsgetn(), and xsputn() methods are
overridden; all the others are inherited as no-ops from std::streambuf.
*/
class rawfilebuf : public sam::streambuf {
public:
  /// Construct a closed buffer
  rawfilebuf() : fd_(-1), owned_(false) { }

  /// @brief Destroy this buffer object, optionally closing the underlying
  /// file descriptor
  ~rawfilebuf() { if (owned_)  close_nothrow(); }

  /// Open a file that will be closed when this buffer is destroyed
  rawfilebuf* open(const char* fname, std::ios_base::openmode mode,
		   int perm = 0664);

  /// Open a file that will be closed when this buffer is destroyed
  rawfilebuf* open(const char* fname, int flags, int perm = 0664);

  /// Associate an open file descriptor that will be closed when this buffer is destroyed
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
  // FIXME Virtual functions that I maybe should be overriding:
  // imbue (default does nothing, ok)
  // overflow (for output; default does nothing, ok)
  // pbackfail (default just fails, ok)
  // setbuf (default does nothing, ok)
  // sync (default does nothing, ok)
  // uflow (default calls underflow(), ok if that returns EOF)
  // underflow (default just returns EOF, ok)


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
