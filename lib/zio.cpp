#include "lib/zio.h"

#include <string>
#include <cstring>

#include <zlib.h>

#include "sam/exception.h"
#include "lib/bgzf.h"
#include "lib/utilities.h"

using std::string;

namespace sam {

inline size_t min(size_t a, size_t b) { return (a < b)? a : b; }

bgzfbuf::bgzfbuf(std::streambuf* sbuf0)
  : sbuf(sbuf0), at_eof(false), buffer(65536), cdata(65536 + BGZF::hsize) {
}

bgzfbuf::~bgzfbuf() {
}

std::streamsize bgzfbuf::sbuf_sgetn(char* buffer, std::streamsize length) {
  std::streamsize n = (at_eof)? 0 : sbuf->sgetn(buffer, length);
  if (n == 0)  at_eof = true;
  return n;
}

string zlib_error_message(const char* function, const z_stream& z) {
  string s = "zlib::";
  s += function;
  s += "() failed";
  if (z.msg)
    s += ": ", s += z.msg;
  return s;
}

// Decompress the specified data into  buffer, discarding whatever may have
// been there previously.
size_t bgzfbuf::buffer_inflate(char* data, size_t length) {
  z_stream z;
  z.zalloc = Z_NULL;
  z.zfree  = Z_NULL;

  z.next_in  = reinterpret_cast<unsigned char*>(data);
  z.avail_in = length;
  if (inflateInit2(&z, -15) != Z_OK)
    throw bad_format(zlib_error_message("inflateInit2", z));

  buffer.clear();
  z.next_out  = reinterpret_cast<unsigned char*>(buffer.begin);
  z.avail_out = buffer.capacity();
  if (inflate(&z, Z_FINISH) != Z_STREAM_END) {
    string message = zlib_error_message("inflate", z);
    (void) inflateEnd(&z);
    throw bad_format(message);
  }

  buffer.end += z.total_out;

  if (inflateEnd(&z) != Z_OK)
    throw bad_format(zlib_error_message("inflateEnd", z));

  return z.total_in;
}

// Read from the streambuf and decompress to fill  buffer, which is assumed
// to be previously empty.  Returns true if the buffer is nonempty afterwards.
// As this is a BGZF file, in fact exactly one BGZF block is decompressed.
bool bgzfbuf::underflow() {
  // Read from the streambuf if necessary to get a complete BGZF header.
  if (cdata.size() < BGZF::hsize) {
    cdata.flush();
    cdata.end += sbuf_sgetn(cdata.end, BGZF::hsize - cdata.size());

    // If there's still no data, we're cleanly at EOF.
    if (cdata.size() == 0)  return false;
  }

  if (BGZF::is_bgzf_header(cdata.begin, cdata.size())) {
    size_t blockonly_length = BGZF::block_size(cdata.begin) - BGZF::hsize;
    cdata.begin += BGZF::hsize;
    cdata.flush();
    // Read the BGZF block and also the following block's header.
    cdata.end += sbuf_sgetn(cdata.end,
			    blockonly_length - cdata.size() + BGZF::hsize);
    if (cdata.size() < blockonly_length)
      throw bad_format(make_string()
	  << "Truncated BGZF block (expected " << blockonly_length
	  << " bytes after header; got " << cdata.size() << ")");

    cdata.begin += buffer_inflate(cdata.begin, blockonly_length);
    return true;
  }
  else
    throw bad_format("Invalid BGZF block header");
}

size_t bgzfbuf::sgetn(char* dest, size_t desired_length) {
  // TODO  Ideally this would unpack just enough of the stream to fill DEST
  // and decompress directly into the destination buffer, thus saving a copy.
  // However, it is easier (and perhaps faster) to decompress an entire block
  // at once.

  size_t length = 0;

  while (true) {
    size_t copy_length = min(desired_length, buffer.end - buffer.begin);
    memcpy(dest, buffer.begin, copy_length);
    buffer.begin += copy_length;
    dest += copy_length;
    length += copy_length;
    desired_length -= copy_length;

    if (desired_length == 0)
      break;
    else if (! underflow())
      break;
  }

  return length;
}

#if 0
void bgzfbuf::sputn(const char* buffer, size_t length) {
}
#endif

} // namespace sam
