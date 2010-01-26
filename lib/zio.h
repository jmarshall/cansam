#ifndef CANSAM_BGZFBUF_H
#define CANSAM_BGZFBUF_H

#include <streambuf>

#include "lib/buffer.h"

namespace sam {

class bgzfbuf {
public:
  bgzfbuf(std::streambuf* sbuf0);
  ~bgzfbuf();

  // Reads up to LENGTH decompressed bytes into BUFFER, returning how many
  // bytes were in fact read.
  size_t sgetn(char* buffer, size_t length);

  void sputn(const char* buffer, size_t length);

private:
  std::streamsize sbuf_sgetn(char* buffer, std::streamsize length);
  bool underflow();
  size_t buffer_inflate(char* data, size_t length);

  std::streambuf* sbuf;
  bool at_eof;

  read_buffer buffer;
  read_buffer cdata;
};

} // namespace sam

#endif
