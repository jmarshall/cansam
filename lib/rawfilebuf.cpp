#include <sys/ioctl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "sam/iobuffer.h"

namespace sam {

// FIXME This paranoia is hopefully unwarranted, and we can just use
// the EOF-returning defaults supplied by std::streambuf.
std::streambuf::int_type rawfilebuf::uflow() {
#ifdef BE_NICE_AND_EFFICIENT
  throw std::ios::failure("rawfilebuf::uflow() invoked");
#else
  char c;
  std::streamsize nread = xsgetn(&c, 1); 
  return (nread == 1)? int_type(c) : EOF;
#endif
}

std::streambuf::int_type rawfilebuf::underflow() {
  throw std::ios::failure("rawfilebuf::underflow() invoked");
}

std::streamsize rawfilebuf::xsgetn(char* s, std::streamsize n) {
  ssize_t nread;
  while ((nread = ::read(fd_, s, n)) < 0)
    if (errno != EINTR)
      throw std::ios::failure("rawfilebuf::xsgetn(): error reading the file");

  return nread;
}

std::streamsize rawfilebuf::xsputn(const char* s, std::streamsize n) {
  std::streamsize total = 0;

  while (n > 0) {
    ssize_t nwritten = ::write(fd_, s, n);
    if (nwritten < 0) {
      if (errno == EINTR)
	nwritten = 0;
      else
	throw std::ios::failure("rawfilebuf::xsputn(): error writing the file");
    }

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

} // namespace sam
