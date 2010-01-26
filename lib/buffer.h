#ifndef CANSAM_BUFFER_H
#define CANSAM_BUFFER_H

#include <cstring>
#include <boost/scoped_array.hpp>

namespace sam {

class read_buffer {
public:
  // Allocate a buffer of the specified capacity
  read_buffer(size_t sz) : array(new char[sz]), capacity_(sz) { clear(); }

  ~read_buffer() { }

  char* begin;
  char* end;

  // Returns the total capacity of the buffer
  size_t capacity() const { return capacity_; }

  // Returns the number of characters currently unread in the buffer
  size_t size() const { return end - begin; }

  // Empty the buffer, updating begin/end to point to its start
  void clear() { begin = end = array.get(); }

  // Move any unread characters to the start of the buffer, updating
  // begin/end accordingly
  void flush() {
    memmove(array.get(), begin, end - begin);
    end = &array[end - begin];
    begin = &array[0];
  }

private:
  boost::scoped_array<char> array;
  size_t capacity_;
};

} // namespace sam

#endif
