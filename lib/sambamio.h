/*  sambamio.h -- SAM/BAM input/output formatting.

    Copyright (C) 2010 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
 3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
    nor the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND ITS CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */

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

  // FIXME isamstream or samstream_base?  is it an isamstream in the ctor already? yes, but are we happy with that?
  static sambamio* new_in(isamstream&);
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
  static std::streamsize
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
