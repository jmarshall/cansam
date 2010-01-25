#ifndef SAMBAMIO_H
#define SAMBAMIO_H

#include "sam/stream.h"

namespace sam {

class collection;
class alignment;

class samstream_base::sambamio {
public:
  class eof_exception { };

  static sambamio* new_in(std::streambuf* sbuf);
  static sambamio* new_out(std::streambuf* sbuf, openmode);

  virtual ~sambamio() { }

  virtual bool get(isamstream&, collection&) = 0;

  // Returns true when an alignment has been successfully read,
  // false at EOF, or throws an exception on formatting or I/O errors.
  virtual bool get(isamstream&, alignment&) = 0;

  virtual void put(osamstream&, const alignment&) = 0;
};

} // namespace sam

#endif
