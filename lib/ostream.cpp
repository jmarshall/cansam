/*  ostream.cpp -- Various output stream inserters.

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

#include <ostream>

#include "sam/alignment.h"
#include "sam/header.h"
#include "lib/utilities.h"  // FIXME NUKE-ME want make_string for dump_on()

using std::string;

namespace {

/* The get_buffer() function uses a std::ios_base iword/pword pair to attach
a buffer to each stream used, for temporary use by the inserters below.
Possible states for this pair, and their implications for get_buffer(), are:

  0, NULL  Initial values; need to allocate buffer and register callback
 -1, NULL  Values after copyfmt(); need to allocate buffer
+ve, ptr   Normal state; may need to reallocate buffer.  */

void callback(std::ios_base::event event, std::ios_base& stream, int index) {
  switch (event) {
  case std::ios_base::erase_event:
    delete [] static_cast<char*>(stream.pword(index));
    break;

  case std::ios_base::copyfmt_event:
    // We can't share the source stream's buffer and would rather not allocate
    // one here (because it is tedious to report failure from a callback), so
    // instead we leave a state that ensures get_buffer() needs to allocate.
    stream.iword(index) = -1;
    stream.pword(index) = NULL;
    break;

  case std::ios_base::imbue_event:
    break;
  }
}

// Returns a buffer of at least the specified capacity local to the stream.
char* get_buffer(std::ios_base& stream, int capacity) {
  static const int index = std::ios_base::xalloc();

  int oldcapacity = stream.iword(index);
  if (capacity > oldcapacity) {
    capacity += 200;  // In case subsequent needs are only a little larger
    char* newbuffer = new char[capacity];

    void*& ptr = stream.pword(index);
    char* oldbuffer = static_cast<char*>(ptr);
    ptr = newbuffer;
    stream.iword(index) = capacity;
    delete [] oldbuffer;

    if (oldcapacity == 0)
      stream.register_callback(callback, index);
  }

  return static_cast<char*>(stream.pword(index));
}

} // anonymous namespace

namespace sam {

std::ostream& operator<< (std::ostream& out, const header& header) {
  char* buffer = get_buffer(out, header.sam_length() + 1);
  *format_sam(buffer, header) = '\0';
  return out << buffer;
}

std::ostream& operator<< (std::ostream& out, const collection& headers) {
  for (collection::const_iterator it = headers.begin();
       it != headers.end(); ++it)
    out << *it << '\n';

  if (out.flags() & std::ios::showpoint) {
    int i = 0;
    out << "Reflist:";
    for (collection::const_ref_iterator it = headers.ref_begin();
	 it != headers.ref_end(); ++it, ++i)
      out << " " << i << "->" << it->name();
    if (! headers.refseqs_in_headers)  out << "  (owned)";
    out << "\nRefmap:";
    for (std::map<std::string, refsequence*>::const_iterator
	 it = headers.refnames.begin(); it != headers.refnames.end(); ++it)
      out << " " << it->first << "->" << it->second->name();
    out << "\n";
  }

  return out;
}

std::ostream& operator<< (std::ostream& out, const header::tagfield& field) {
#ifdef STILL_TAB_DELIMITED
  const char* limit = header::tagfield::next(field.tag_);
  return out << string(field.tag_, limit - field.tag_);
#else
  return out << field.tag_;
#endif
}

std::ostream& operator<< (std::ostream& out, header::const_iterator it) {
  return out << static_cast<const void*>(it.ptr);
}

std::ostream& operator<< (std::ostream& out, const alignment& aln) {
  char* buffer = get_buffer(out, aln.sam_length() + 1);
  *format_sam(buffer, aln, out) = '\0';
  return out << buffer;
}

std::ostream& operator<< (std::ostream& out, const alignment::tagfield& aux) {
  char* buffer = get_buffer(out, aux.sam_length() + 1);
  *format_sam(buffer, aux) = '\0';
  return out << buffer;
}

std::ostream& operator<< (std::ostream& out, alignment::const_iterator it) {
  return out << static_cast<const void*>(it.ptr);
}

// FIXME nuke me or otherwise hide me!
void alignment::dump_on(std::ostream& out, const_iterator marker) const {
  make_string text;
  const char* s = p->data();
  const char* limit = &s[p->capacity()];
  while (s < limit) {
    if (s == qname_c_str())  text << "]NAME:[";
    if (s == p->cigar_data())  text << "]CIG:[";
    if (s == seq_raw_data())  text << "]SEQ:[";
    if (s == qual_raw_data())  text << "]QUAL:[";
    if (s == begin().ptr)  text << "]AUXEN:[";

    if (s == marker.ptr)  text << " [" << *s++ << "] ";
    else  text << *s++;
  }

  out << "Capacity:" << p->capacity() << ", cindex:" << p->h.cindex
      << ", data:[" << string(text) << "]\n";
}

} // namespace sam
