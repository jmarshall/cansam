#include "sam/streambuf.h"

#include <stdexcept>

#include <sys/ioctl.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "sam/exception.h"

namespace sam {

rawfilebuf*
rawfilebuf::open(const char* fname, std::ios_base::openmode mode, int perm) {
  using std::ios;

  // FIXME Either nuke or prettify table
  //
  //    out		w	O_WRONLY|O_CREAT|O_TRUNC
  //    out trunc	w	O_WRONLY|O_CREAT|O_TRUNC
  //    out       app	a	O_WRONLY|O_CREAT        |O_APPEND
  // in			r	O_RDONLY
  // in out		r+	O_RDWR
  // in out trunc	w+	O_RDWR  |O_CREAT|O_TRUNC
  // in out       app	a+	O_RDWR  |O_CREAT|       |O_APPEND    (not std)
  //
  // (Actually will be std, as will in|app, app, as per DR 596.)

  int flags = 0;
  if (mode & ios::in) {
    flags = (mode & ios::out)? O_RDWR : O_RDONLY;
  }
  else if (mode & ios::out)  {
    flags = O_WRONLY | O_CREAT;
    if (! (mode & ios::app))  flags |= O_TRUNC;
  }

  if (mode & ios::trunc)   flags |= O_CREAT | O_TRUNC;
  if (mode & ios::app)     flags |= O_CREAT | O_APPEND;
#ifdef O_BINARY
  if (mode & ios::binary)  flags |= O_BINARY;
#endif
#ifdef O_TEXT
  // FIXME  Actually probably not a good idea, at least on Cygwin
  if (! (mode & ios::binary))  flags |= O_TEXT;
#endif

  if (! open(fname, flags, perm))
    return NULL;

  if (mode & ios::ate) {
    if (::lseek(fd_, 0, SEEK_END) <= 0) {
      int saved_errno = errno;
      close_nothrow();
      errno = saved_errno;
      return NULL;
    }
  }

  return this;
}

rawfilebuf* rawfilebuf::open(const char* fname, int flags, int perm) {
  if (is_open())  return NULL;

  do fd_ = ::open(fname, flags, perm); while (fd_ < 0 && errno == EINTR);
  if (fd_ < 0)
    return NULL;

  owned_ = true;
  return this;
}

// Close the underlying file descriptor, failing (i.e., returning negative)
// on hard errors rather than throwing an exception.
int rawfilebuf::close_nothrow() {
  if (! is_open())  return 0;

  int ret;
  do ret = ::close(fd_); while (ret < 0 && errno == EINTR);

  fd_ = -1;
  return ret;
}

void rawfilebuf::close() {
  if (close_nothrow() < 0)
    throw sam::system_error("close() failed", errno);
}

std::streamsize rawfilebuf::xsgetn(char* s, std::streamsize n) {
  ssize_t nread;
  do nread = ::read(fd_, s, n); while (nread < 0 && errno == EINTR);
  if (nread < 0)
    throw sam::system_error("read() failed", errno);

  return nread;
}

std::streamsize rawfilebuf::xsputn(const char* s, std::streamsize n) {
  std::streamsize total = 0;

  while (n > 0) {
    ssize_t nwritten;
    do nwritten = ::write(fd_, s, n); while (nwritten < 0 && errno == EINTR);
    if (nwritten < 0)
      throw sam::system_error("write() failed", errno);

    total += nwritten;
    s += nwritten;
    n -= nwritten;
  }

  return total;
}

std::streamsize rawfilebuf::showmanyc() {
  off_t pos = ::lseek(fd_, 0, SEEK_CUR);
  if (pos >= 0) {
    struct stat st;
    if (::fstat(fd_, &st) == 0)
      return st.st_size - pos;
  }
  else if (errno == ESPIPE) {
    // We're not seekable, so perhaps we're STREAMS-based or otherwise
    // responsive to this ioctl.
#ifdef FIONREAD
    int n = 0;
    if (::ioctl(fd_, FIONREAD, &n) == 0 && n >= 0)
      return n;
#endif
  }

  return 0;
}

std::streampos
rawfilebuf::seekoff(std::streamoff off, std::ios_base::seekdir way,
		    std::ios_base::openmode) {
  int whence = (way == std::ios::beg)? SEEK_SET :
	       (way == std::ios::cur)? SEEK_CUR : SEEK_END;

  return ::lseek(fd_, off, whence);
}

std::streampos
rawfilebuf::seekpos(std::streampos pos, std::ios_base::openmode) {
  return ::lseek(fd_, pos, SEEK_SET);
}

std::streambuf::int_type rawfilebuf::uflow() {
  throw std::logic_error("rawfilebuf::uflow() invoked");
}

std::streambuf::int_type rawfilebuf::underflow() {
  throw std::logic_error("rawfilebuf::underflow() invoked");
}

std::streambuf::int_type rawfilebuf::overflow(int_type) {
  throw std::logic_error("rawfilebuf::overflow() invoked");
}

} // namespace sam
