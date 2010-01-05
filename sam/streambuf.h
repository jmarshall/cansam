/** @file sam/streambuf.h
    @brief Low-level classes for input/output

Provides low-level input/output classes derived from @c std::streambuf.
Most code will not need to use these classes directly.  */

#ifndef CANSAM_STREAMBUF_H
#define CANSAM_STREAMBUF_H

#include <streambuf>

namespace sam {

/** @class sam::streambuf sam/streambuf.h
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

/** @class sam::rawfilebuf sam/streambuf.h
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
