#ifndef SAMBAMIO_H
#define SAMBAMIO_H

#include <vector>

#include "sam/stream.h"

namespace sam {

class alignment;
class char_buffer;
class collection;

class sambamio {
public:
  class eof_exception { };

  static sambamio* new_in(std::streambuf*);
  static sambamio* new_out(std::ios::openmode);

  virtual ~sambamio() { }

  virtual bool get(isamstream&, collection&) = 0;

  // Returns true when an alignment has been successfully read,
  // false at EOF, or throws an exception on formatting or I/O errors.
  virtual bool get(isamstream&, alignment&) = 0;

  virtual void put(osamstream&, const collection&) = 0;
  virtual void put(osamstream&, const alignment&) = 0;
  virtual void flush(osamstream&) = 0;

protected:
  sambamio() : header_cindex(0) { }

  // FIXME blah blah
  std::streamsize
  rdbuf_sgetn(isamstream& stream, char* buffer, std::streamsize length) {
    // FIXME can we get rid of this by always writing it out in the caller...?
    if (stream.eof())  return 0;

    std::streamsize n = stream.rdbuf()->sgetn(buffer, length);
    if (n == 0) {
      if (stream.setstate_wouldthrow(std::ios::eofbit))  throw eof_exception();
    }
    return n;
  }

  // Initialise BUFFER for use with peek() & getline(), optionally prefilling
  // with TEXT/TEXTSIZE.  These functions will call xsgetn() to refill BUFFER
  // as needed.
  void prepare_line_buffer(char_buffer& buffer,
			   const char* text = NULL, size_t textsize = 0);

  int peek(char_buffer&, isamstream&);
  int getline(char_buffer&, isamstream&, std::vector<char*>&);

  // For use by peek() & getline().
  virtual size_t xsgetn(isamstream&, char*, size_t) = 0;

  // For use by get(collection&).
  void set_cindex(collection&);

  // Cached copy of the header's cindex, for use by get(alignment&).
  int header_cindex;
};

// Bitmask flags for use with the private collection::push_back().
enum { add_header = 1, add_refseq = 2, add_refname = 4 };

} // namespace sam

#endif
