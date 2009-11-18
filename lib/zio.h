#ifndef CANSAM_FLATE_H
#define CANSAM_FLATE_H

#include <streambuf>

namespace sam {
namespace internal {

class flatebuf {
public:
  flatebuf(std::streambuf* sbuf) : sbuf_(sbuf) { }
  ~flatebuf() { }

  // Reads up to LENGTH decompressed bytes into BUFFER, returning how many
  // bytes were in fact read.
  size_t read_uncompressed(char* buffer, size_t length);

  void write_compressed(const char* buffer, size_t length);

private:
  std::streambuf* sbuf_;
};

} // namespace internal
} // namespace sam

#endif
