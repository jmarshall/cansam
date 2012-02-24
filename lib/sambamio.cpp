/*  sambamio.cpp -- SAM/BAM input/output formatting.

    Copyright (C) 2010-2012 Genome Research Ltd.

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

#include "lib/sambamio.h"

#include <streambuf>
#include <string>
#include <utility>
#include <vector>
#include <cstddef>
#include <cstring>

#include <iostream> // FIXME NUKE-ME

#include <zlib.h>

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"
#include "cansam/sam/stream.h"
#include "cansam/exception.h"
#include "lib/utilities.h"
#include "lib/wire.h"

#ifndef EOF
#define EOF  (std::char_traits<char>::eof())
#endif

using std::string;

namespace sam {

namespace BGZF {
/* BAM files consist of BGZF blocks, which are RFC 1952 GZIP members with
a "BC" extra subfield in their headers, typically as follows:

    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17   Offset from start
   1f 8b .. *4 .. .. .. .. .. .. 06 00 42 43 02 00 .. ..   BGZF signature
   magic    ^bit 2 is FEXTRA flag      ^B ^C

These BGZF utilities provide what we need for identifying GZIP and BGZF headers
and extracting the interesting field from the "BC" subfield.  */

// Size, in bytes, of a BGZF block header
enum { hsize = 18, tsize = 8, full_block_size = 65536,
       payload_max_size = full_block_size - hsize - tsize };

enum { uncompressed_max_size = 65536 };

// Returns whether the specified memory block starts with a GZIP member header
inline bool is_gzip_header(const char* s, int length) {
  return length >= 2 && s[0] == '\x1f' && s[1] == '\x8b';
}

// Returns whether the specified memory block starts with a valid BGZF header
inline bool is_bgzf_header(const char* s, int length) {
  return length >= 18 && s[0] == '\x1f' && s[1] == '\x8b' &&
	 (s[3] & 4) && memcmp(&s[10], "\6\0\x42\x43\2\0", 6) == 0;
}

// For a valid BGZF header, returns the block_size field
inline int block_size(const char* s) {
  return convert::uint16(&s[16]) + 1;
}

// Write a BGZF block header, and return the number of bytes written
inline int write_bgzf_header(char* s, int block_size) {
  static const char boilerplate[] =
    { '\x1f', '\x8b', 8, '\x04', 0,0,0,0, 0, '\xff', 6,0, '\x42','\x43', 2,0 };

  memcpy(s, boilerplate, sizeof boilerplate);
  convert::set_bam_uint16(&s[16], block_size - 1);
  return hsize;
}

// Write a BGZF block trailer, and return the number of bytes written
inline int write_bgzf_trailer(char* s, uint32_t crc, int uncompressed_size) {
  convert::set_bam_uint32(s, crc);
  convert::set_bam_uint32(&s[4], uncompressed_size);
  return tsize;
}

} // namespace BGZF

class char_buffer {
public:
  // Allocate a buffer with the specified capacity
  char_buffer(size_t sz) : array(new char[sz]), capacity(sz) { clear(); }

  ~char_buffer() { delete [] array; }

  char* begin;
  char* end;

  // Returns the number of characters currently unread in the buffer
  size_t size() const { return end - begin; }

  // Returns the number of character positions available beyond  end.
  size_t available() const { return &array[capacity] - end; }

  // Empty the buffer, updating begin/end to point to its start
  void clear() { begin = end = array; }

  // Move any unread characters to the start of the buffer, updating
  // begin/end accordingly
  void flush() {
    memmove(array, begin, end - begin);
    end = &array[end - begin];
    begin = array;
  }

  // Ensure the buffer has at least the specified caacity
  void reserve(size_t sz) { if (capacity < sz)  reserve_(sz); }

  // Move any unread characters to the start of the buffer and ensure there is
  // reasonable space available beyond  end, possibly by enlarging the buffer.
  // Updates begin/end and the specified pointers accordingly.
  void flush_make_space(char*& ptr, std::vector<char*>& ptrvec);

private:
  void reserve_(size_t sz);

  char* array;
  size_t capacity;
};

void char_buffer::reserve_(size_t sz) {
  char* oldarray = array;
  char* oldbegin = begin;

  size_t newcapacity = sz + 32768;
  char* newarray = new char[newcapacity];
  memcpy(newarray, begin, end - begin);
  array = newarray;
  capacity = newcapacity;

  begin = array;
  end = &begin[end - oldbegin];

  delete [] oldarray;
}

void char_buffer::flush_make_space(char*& ptr, std::vector<char*>& ptrvec) {
  char* oldarray = NULL;
  char* oldbegin = begin;

  if (begin > array) {
    memmove(array, begin, end - begin);
  }
  else if (available() < 1024 + capacity / 4) {
    size_t newcapacity = capacity * 2;
    char* newarray = new char[newcapacity];
    oldarray = array;
    memcpy(newarray, begin, end - begin);
    array = newarray;
    capacity = newcapacity;
  }
  else
    return;  // Don't bother rewiring the pointers when nothing has changed.

  begin = array;
  end = &begin[end - oldbegin];
  ptr = &begin[ptr - oldbegin];

  for (std::vector<char*>::iterator it = ptrvec.begin();
       it != ptrvec.end(); ++it)
    *it = &begin[*it - oldbegin];

  delete [] oldarray;
}


/* Pointers to the read buffer are arranged as follows:

  [---------ABCDEF/GHIJKL/MNO
  ^         ^                ^
  buffer    begin            end */


void sambamio::prepare_line_buffer(char_buffer& b,
				   const char* text, size_t textsize) {
  if (text) {
    memcpy(b.end, text, textsize);
    b.end += textsize;
  }

  *b.end = '\n';  // Set up us the sentinel.
}

int sambamio::peek(char_buffer& b, isamstream& stream) {
  if (b.begin < b.end)
    return *b.begin;

  if (! stream.eof()) {
    b.flush();
    // Read more characters, leaving one position spare for the sentinel.
    b.end += xsgetn(stream, b.end, b.available() - 1);
    *b.end = '\n';
  }

  return (b.begin < b.end)? *b.begin : EOF;
}

/* Reads a newline-terminated line of tab-delimited text into  fields,
and returns the number of fields present (or 0 at EOF).
*/
int sambamio::getline(char_buffer& b, isamstream& stream,
		      std::vector<char*>& fields) {
  fields.clear();
  fields.push_back(b.begin);

  char* s = b.begin;
  while (true)
    if (*s == '\t') {
      *s++ = '\0';
      fields.push_back(s);
    }
    else if (*s == '\n') {
      if (s < b.end) {
	// This is a real \n, so a properly-terminated line has been read.
	if (s > b.begin && s[-1] == '\r')  s[-1] = '\0', fields.push_back(s++);
	else  *s++ = '\0', fields.push_back(s);

	break;
      }
      else if (stream.eof()) {
	// This is the sentinel and no further characters are forthcoming.
	// If any characters have been read, they constitute a final line
	// (which is unterminated); otherwise we are properly at EOF.
	if (s > b.begin) {
	  // Move the sentinel one character later, to make room for a \0.
	  b.end++;
	  if (b.available() == 0)  b.flush_make_space(s, fields);
	  *b.end = '\n';

	  *s++ = '\0';
	  fields.push_back(s);
	}

	break;
      }
      else {
	// This is the sentinel, and there are more characters to be read.
	b.flush_make_space(s, fields);

	// Read more characters, leaving one position spare for the sentinel.
	b.end += xsgetn(stream, b.end, b.available() - 1);
	*b.end = '\n';
      }
    }
    else
      s++;

  b.begin = s;
  return fields.size() - 1;
}

// Called from get(collection&); allocates a new cindex, as we're effectively
// assigning a new collection that's unrelated to whatever was previously
// in  headers.  Caches the new cindex for use in get(alignment&).
void sambamio::set_cindex(collection& headers) {
  headers.reallocate_cindex();
  header_cindex = headers.cindex;
}

inline size_t min(size_t a, size_t b) { return (a < b)? a : b; }

/* A sam::alignment object contains only a pointer to a variable-sized memory
block containing all the alignment data represented mostly in the same way as
it is represented in an (uncompressed) BAM file.

  +---------------+---------+--...--+-...-+-...--+--...--+------...------+
  | capacity, etc | bamcore | cigar | seq | qual | rname | aux fields... |
  +---------------+---------+--...--+-...-+-...--+--...--+------...------+

The memory block is allocated as a suitably-aligned appropriately-sized char
array and then interpreted as a sam::alignment::block.  This block contains
some housekeeping data followed by the alignment data as it appears in a BAM
file, except that the variable-length fields have been rearranged so that the
packed cigar and seq fields are first and thus naturally aligned.  */


// Binary BAM files
// ================

class bamio : public sambamio {
public:
  bamio(const char* text, std::streamsize textsize);
  bamio(bool compression);
  virtual ~bamio();

  virtual bool get(isamstream&, collection&);
  virtual bool get(isamstream&, alignment&);
  virtual void put(osamstream&, const collection&);
  virtual void put(osamstream&, const alignment&);
  virtual void flush(osamstream&);

protected:
  virtual size_t xsgetn(isamstream&, char*, size_t);

private:
  size_t read(isamstream&, void*, size_t);
  int32_t read_int32(isamstream&);
  void read_refinfo(isamstream& stream, string& name, coord_t& length);

  bool underflow(isamstream&);
  void fill_cdata(isamstream&, size_t);
  size_t inflate_into_buffer(char*, size_t);

  size_t deflate_onto_cdata(osamstream&, char*, size_t);

  static string zlib_message(const char* function, const z_stream& z);

  char_buffer buffer;
  char_buffer cdata;

  int compression_level;      // Used in deflate_onto_cdata()
  size_t header_text_length;  // Used in xsgetn()

  z_stream zinflate, zdeflate;
  bool zinflate_active, zdeflate_active;

  typedef unsigned char uchar;  // For casting to the pointers exposed by zlib

  // Hmmm...
  void flush_buffer(osamstream&);
};

// Called when blah blah blah
void bamio::flush_buffer(osamstream& stream) {
  buffer.begin += deflate_onto_cdata(stream, buffer.begin, min(buffer.size(), BGZF::uncompressed_max_size));
  buffer.flush();
}

bamio::bamio(const char* text, std::streamsize textsize)
  : buffer(65536), cdata(65536),
    zinflate_active(false), zdeflate_active(false) {
  memcpy(cdata.end, text, textsize);
  cdata.end += textsize;
}

bamio::bamio(bool compression)
  : buffer(65536), cdata(65536),
    compression_level(compression? Z_DEFAULT_COMPRESSION : Z_NO_COMPRESSION),
    zinflate_active(false), zdeflate_active(false) {
}

string bamio::zlib_message(const char* function, const z_stream& z) {
  make_string s;
  s << "zlib::" << function << "() failed";
  if (z.msg)  s << ": " << z.msg;
  return s;
}

bamio::~bamio() {
  // In general, destructors should not throw exceptions.  Thus Z_DATA_ERROR
  // is ignored below, as in this case the problem will already have been
  // reported by an earlier zlib function and caused an exception to be
  // thrown, so the present code will likely have been reached during the
  // resulting stack unwinding -- hence must not itself throw.
  //
  // However Z_STREAM_ERRORs reported by these functions represent bugs in
  // the calling program; since this Can't Happen, we'd like to hear about
  // it, either by exception or via terminate().

  if (zinflate_active) {
    if (inflateEnd(&zinflate) != Z_OK)
      throw std::logic_error(zlib_message("inflateEnd", zinflate));
  }

  if (zdeflate_active) {
    int status = deflateEnd(&zdeflate);
    if (status != Z_OK && status != Z_DATA_ERROR)
      throw std::logic_error(zlib_message("deflateEnd", zdeflate));
  }
}

// Fill  cdata  by reading from the streambuf.  Reads at least  desired_size
// bytes (unless EOF occurs first), but potentially as much as will fit into
// the buffer, which is probably several BGZF blocks.  This reduces the number
// of read(2) system calls for sequential access, and (with modern disk block
// sizes) shouldn't mean too much data is discarded by any subsequent seek().
void bamio::fill_cdata(isamstream& stream, size_t desired_size) {
  cdata.flush();

  if (! stream.eof())
    do {
      std::streamsize n = rdbuf_sgetn(stream, cdata.end, cdata.available());
      if (n == 0)  break;
      cdata.end += n;
    } while (cdata.size() < desired_size);
}

// Compress the specified data into  cdata, appending to whatever is already
// in  cdata  and writing that to the stream if necessary to make space.
// Returns the number of input bytes actually compressed into the output block.
size_t
bamio::deflate_onto_cdata(osamstream& stream, char* data, size_t length) {
  while (true) {
    zdeflate.next_in = reinterpret_cast<uchar*>(data);
    zdeflate.avail_in = length;

    if (zdeflate_active) {
      if (deflateReset(&zdeflate) != Z_OK)
	throw std::logic_error(zlib_message("deflateReset", zdeflate));
    }
    else {
      zdeflate.zalloc = Z_NULL;
      zdeflate.zfree  = Z_NULL;
      if (deflateInit2(&zdeflate, compression_level, Z_DEFLATED,
		       -15, MAX_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK)
	throw std::logic_error(zlib_message("deflateInit2", zdeflate));
      zdeflate_active = true;
    }

    // At this point, and in general, the invariant is that  cdata.[begin,end)
    // contains buffered previously-deflated blocks.  We compress the new block
    // into the available space after  cdata.end.

    if (cdata.available() < 10 * 1024) {
      // Ensure that there is space available for the block header
      // and a reasonable stab at compression.
      while (cdata.size() > 0)
	cdata.begin += stream.rdbuf()->sputn(cdata.begin, cdata.size());
      cdata.flush();
    }

    zdeflate.next_out = reinterpret_cast<uchar*>(cdata.end + BGZF::hsize);
    zdeflate.avail_out = min(cdata.available() - BGZF::hsize - BGZF::tsize,
			     BGZF::payload_max_size);

    int status = deflate(&zdeflate, Z_FINISH);
    if (status == Z_OK) {
      // There was not enough output space, either because there really isn't
      // much space available in  cdata  or (conceivably) because the data
      // cannot be compressed within the block payload size.  We flush previous
      // blocks from  cdata  and continue deflate()ing -- the capacity of  cdata
      // exceeds the BGZF block size, so after flushing these two cases can be
      // distinguished.

      while (cdata.size() > 0)
	cdata.begin += stream.rdbuf()->sputn(cdata.begin, cdata.size());

      // Temporarily extend  cdata.[begin,end)  to include the partially-
      // written block, so that it will be maintained across the flush().
      cdata.end = reinterpret_cast<char*>(zdeflate.next_out);
      cdata.flush();

      zdeflate.next_out = reinterpret_cast<uchar*>(cdata.end);
      zdeflate.avail_out = BGZF::payload_max_size - zdeflate.total_out;

      // Re-establish the general  cdata.[begin,end)  invariant.
      cdata.end = cdata.begin;

      if (zdeflate.avail_out > 0)
	status = deflate(&zdeflate, Z_FINISH);
    }

    // Normally deflation was successfully completed, and we break out of
    // this loop the first time through...
    if (status == Z_STREAM_END)
      break;
    else if (status != Z_OK)
      throw sam::exception(zlib_message("deflate", zdeflate));

    // ...otherwise this is the unusual case of incompressible data that
    // cannot be made to fit within the block payload size.  Try again, but
    // only attempt to compress approximately the amount of data known to fit.
    if (zdeflate.total_in >= 1024)
      length = zdeflate.total_in - 128;
    else
      throw std::logic_error("implausibly incompressible data");
  }

  uint32_t crc =
      crc32(crc32(0, NULL, 0), reinterpret_cast<uchar*>(data), length);

  // Fill in the BGZF header and trailer, and extend  cdata.[begin,end)
  // to encompass the newly buffered block.
  cdata.end += BGZF::write_bgzf_header(cdata.end,
			BGZF::hsize + zdeflate.total_out + BGZF::tsize);
  cdata.end += zdeflate.total_out;
  cdata.end += BGZF::write_bgzf_trailer(cdata.end, crc, length);
  return length;
}

// Decompress the specified data into  buffer, discarding whatever may have
// been there previously.
// FIXME Maybe don't have underflow() skip the header so we can do
// inflateInit() normal not raw mode and get CRC checking
size_t bamio::inflate_into_buffer(char* data, size_t length) {
  zinflate.next_in  = reinterpret_cast<uchar*>(data);
  zinflate.avail_in = length;

  if (zinflate_active) {
    if (inflateReset(&zinflate) != Z_OK)
      throw std::logic_error(zlib_message("inflateReset", zinflate));
  }
  else {
    zinflate.zalloc = Z_NULL;
    zinflate.zfree  = Z_NULL;
    if (inflateInit2(&zinflate, -15) != Z_OK)
      throw std::logic_error(zlib_message("inflateInit2", zinflate));
    zinflate_active = true;
  }

  buffer.clear();
  zinflate.next_out  = reinterpret_cast<uchar*>(buffer.end);
  zinflate.avail_out = buffer.available();
  if (inflate(&zinflate, Z_FINISH) != Z_STREAM_END)
    throw bad_format(zlib_message("inflate", zinflate));

  buffer.end += zinflate.total_out;
  return zinflate.total_in;
}

// Decompress one BGZF block into  buffer,  which is assumed to be previously
// empty, from  cdata,  which is refilled by reading from the streambuf if
// necessary.  Returns true if  buffer  is nonempty afterwards.
bool bamio::underflow(isamstream& stream) {
  if (cdata.size() < BGZF::hsize) {
    fill_cdata(stream, BGZF::hsize);
    // If there's still no data, we're cleanly at EOF.
    if (cdata.size() == 0)  return false;
  }

  if (BGZF::is_bgzf_header(cdata.begin, cdata.size())) {
    // "Block only" measures compressed-payload & trailer, excluding the header.
    size_t blockonly_length = BGZF::block_size(cdata.begin) - BGZF::hsize;
    cdata.begin += BGZF::hsize;
    if (cdata.size() < blockonly_length) {
      fill_cdata(stream, blockonly_length);
      if (cdata.size() < blockonly_length)
	throw bad_format(make_string()
	    << "Truncated BGZF block (expected " << blockonly_length
	    << " bytes after header; got " << cdata.size() << ")");
    }

    cdata.begin += inflate_into_buffer(cdata.begin, blockonly_length);
    cdata.begin += 8; // FIXME skip the GZIP member footer... should check it

    return true;
  }
  else
    throw bad_format("Invalid BGZF block header");
}

size_t bamio::read(isamstream& stream, void* destv, size_t desired_length) {
  char* dest = static_cast<char*>(destv);

  // TODO  Ideally this would unpack just enough of the stream to fill DEST
  // and decompress directly into the destination buffer, thus saving a copy.
  // However, for now it is easier (and perhaps faster) to decompress an
  // entire block at once.

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
    else if (! underflow(stream))
      break;
  }

  return length;
}

int32_t bamio::read_int32(isamstream& stream) {
  int32_t x;
  if (read(stream, &x, sizeof x) < sizeof x)
    throw bad_format("Truncated BAM header");

  convert::set_int32(x);
  return x;
}

void bamio::read_refinfo(isamstream& stream, string& name, coord_t& length) {
  size_t name_length = read_int32(stream);

  // FIXME We can do better than this...
  // ...but might have to invent better abstraction than read()

  char namebuf[300];
  if (read(stream, namebuf, name_length) < name_length)
    throw bad_format("Truncated BAM header (in reference list)");

  name.assign(namebuf, name_length - 1);
  length = read_int32(stream);
}

size_t bamio::xsgetn(isamstream& stream, char* buffer, size_t length) {
  if (header_text_length > 0) {
    size_t n = read(stream, buffer, min(length, header_text_length));
    header_text_length -= n;
    return n;
  }
  else {
    buffer[0] = '\0';
    return 1;
  }
}

bool bamio::get(isamstream& stream, collection& headers) {
  char magic[4];
  if (read(stream, magic, sizeof magic) < sizeof magic)
    throw bad_format("Truncated BAM magic number");
  if (memcmp(magic, "BAM\1", sizeof magic) != 0)
    throw bad_format(make_string()
	<< "Invalid BAM magic number ('"
	<< magic[0] << magic[1] << magic[2] << magic[3] << "')");

  headers.clear();
  headers.refseqs_in_headers = false;
  set_cindex(headers);

  header_text_length = read_int32(stream);

  char_buffer header_text_buffer(32768);
  prepare_line_buffer(header_text_buffer);
  std::vector<char*> fields;
  while (peek(header_text_buffer, stream) == '@') {
    int nfields = getline(header_text_buffer, stream, fields);
    headers.push_back(string(fields[0], fields[nfields] - fields[0] - 1),
		      add_header | add_refname);
//std::clog << *(headers.headers.back()) << '\n';
  }
//std::clog << "bamio::get: headers: " << headers.headers.size() << "; refseqs: " << headers.refseqs.size() << "; refnames: " << headers.refnames.size() << "\n";
  // There are two cases, depending on whether there are any @SQ text headers.
  headers.refseqs_in_headers = ! headers.refnames.empty();

  // FIXME Is refseqs_in_headers set up right if an exception is thrown in
  // the loops below?  Is it even possible?

  if (headers.refseqs_in_headers) {
    string name;
    coord_t length;
    int ref_count = read_int32(stream);
//std::clog << "# " << ref_count << " binary entries\n";
    for (int index = 0; index < ref_count; index++) {
      read_refinfo(stream, name, length);
//std::clog << "@SQ\tSN:" << name << "\tLN:" << length << '\n';
      // FIXME more checking...
      refsequence* rhdr = headers.refnames[name];
      rhdr->index_ = index;
      headers.refseqs.push_back(rhdr);
    }
  }
  else {
    string name;
    coord_t length;
    int ref_count = read_int32(stream);
    for (int index = 0; index < ref_count; index++) {
      read_refinfo(stream, name, length);
      refsequence* refp = new refsequence(name, length, index);
      if (! headers.refnames.insert(make_pair(name, refp)).second) {
	delete refp;
	throw sam::exception(make_string()
	    << "Reference \"" << name << "\" duplicated in BAM reference list");
      }

      headers.refseqs.push_back(refp);
    }
  }

  return true;
}

bool bamio::get(isamstream& stream, alignment& aln) {
  uint32_t rest_length;

  size_t n = read(stream, &rest_length, sizeof rest_length);
  if (n < sizeof rest_length) {
    if (n == 0)  return false;
    else  throw bad_format("Truncated BAM alignment record");
  }

  convert::set_uint32(rest_length);

  int size = sizeof(rest_length) + rest_length;
  if (aln.p->capacity() < size)
    aln.resize_unshare_discard(size);

  n = read(stream, &aln.p->c.rindex, rest_length);
  if (n < rest_length)
    throw bad_format(make_string()
	<< "Truncated BAM alignment record (got " << n
	<< " bytes of an expected remainder of " << rest_length << ")");

  aln.p->h.cindex = header_cindex;
  aln.p->c.rest_length = rest_length;
  convert::set_int32(aln.p->c.rindex);
  convert::set_int32(aln.p->c.zpos);
  convert::set_uint16(aln.p->c.bin);
  convert::set_uint16(aln.p->c.cigar_length);
  convert::set_uint16(aln.p->c.flags);
  convert::set_int32(aln.p->c.read_length);
  convert::set_int32(aln.p->c.mate_rindex);
  convert::set_int32(aln.p->c.mate_zpos);
  convert::set_int32(aln.p->c.isize);
  return true;
}

void bamio::flush(osamstream& stream) {
#if ONE_OF_THESE
  while (buffer.size() > 0)
    flush_buffer(stream);
#else
  while (buffer.size() > 0)
    buffer.begin += deflate_onto_cdata(stream, buffer.begin,
			  min(buffer.size(), BGZF::uncompressed_max_size));
  buffer.clear();
#endif

  while (cdata.size() > 0)
    cdata.begin += stream.rdbuf()->sputn(cdata.begin, cdata.size());
  cdata.clear();
}

void bamio::put(osamstream& stream, const collection& coln) {
  int header_length = 0;
  for (collection::const_iterator it = coln.begin(); it != coln.end(); ++it)
    header_length += it->sam_length() + 1;

  memcpy(buffer.end, "BAM\1", 4);
  buffer.end += 4;
  convert::set_bam_int32(buffer.end, header_length);
  buffer.end += sizeof(int32_t);

  for (collection::const_iterator it = coln.begin(); it != coln.end(); ++it)
    if (buffer.available() >= size_t(it->sam_length() + 1)) {
      buffer.end = format_sam(buffer.end, *it);
      *buffer.end++ = '\n';
    }
    else {
      // FIXME bamio::put(headers) buffering
      throw std::logic_error("bamio::put(headers) buffering not implemented");
    }

  if (buffer.size() >= BGZF::uncompressed_max_size)
    flush_buffer(stream);

  convert::set_bam_int32(buffer.end, coln.ref_size());
  buffer.end += sizeof(int32_t);

  for (collection::const_ref_iterator it = coln.ref_begin();
       it != coln.ref_end(); ++it) {
    const string& name = it->name();
    int namelen = name.length();

    if (buffer.available() >= sizeof(int32_t) + namelen+1 + sizeof(int32_t)) {
      convert::set_bam_int32(buffer.end, namelen+1);
      buffer.end += sizeof(int32_t);
      memcpy(buffer.end, name.data(), namelen);
      buffer.end += namelen;
      *buffer.end++ = '\0';
      convert::set_bam_int32(buffer.end, it->length());
      buffer.end += sizeof(int32_t);
    }
    else {
      // FIXME bamio::put(headers) buffering
      throw std::logic_error("bamio::put(headers) buffering not implemented");
    }
  }

  if (buffer.size() >= BGZF::uncompressed_max_size)
    flush_buffer(stream);
}

void bamio::put(osamstream& stream, const alignment& aln) {
  aln.sync();

  int length = min(aln.p->size(), buffer.available());
  memcpy(buffer.end, aln.p->data(), length);

  // Because  buffer  has a capacity that exceeds the BGZF uncompressed block
  // size by at least sizeof(bamcore), it is guaranteed that all of this data
  // that needs to be converted is indeed in this first memcopied block.
  convert::set_bam32(buffer.end + offsetof(alignment::bamcore, rest_length));
  convert::set_bam32(buffer.end + offsetof(alignment::bamcore, rindex));
  convert::set_bam32(buffer.end + offsetof(alignment::bamcore, zpos));
  convert::set_bam16(buffer.end + offsetof(alignment::bamcore, bin));
  convert::set_bam16(buffer.end + offsetof(alignment::bamcore, cigar_length));
  convert::set_bam16(buffer.end + offsetof(alignment::bamcore, flags));
  convert::set_bam32(buffer.end + offsetof(alignment::bamcore, read_length));
  convert::set_bam32(buffer.end + offsetof(alignment::bamcore, mate_rindex));
  convert::set_bam32(buffer.end + offsetof(alignment::bamcore, mate_zpos));
  convert::set_bam32(buffer.end + offsetof(alignment::bamcore, isize));
  buffer.end += length;

  // Ensure that the whole alignment record has been copied -- if it has
  // not, it must be because the buffer has been filled.
  int copied = length;
  while (buffer.size() >= BGZF::uncompressed_max_size) {
#if ONE_OF_THESE
    flush_buffer(stream);
#else
    buffer.begin += deflate_onto_cdata(stream, buffer.begin,
				       BGZF::uncompressed_max_size);
    buffer.flush();
#endif

    if (copied < aln.p->size()) {
      length = min(aln.p->size() - copied, buffer.available());
      memcpy(buffer.end, aln.p->data() + copied, length);
      buffer.end += length;
      copied += length;
    }
  }
}


// Text SAM files
// ==============

class samio : public sambamio {
public:
  samio();
  samio(const char* text, std::streamsize textsize);
  virtual ~samio();

  virtual bool get(isamstream&, collection&);
  virtual bool get(isamstream&, alignment&);
  virtual void put(osamstream&, const collection&);
  virtual void put(osamstream&, const alignment&);
  virtual void flush(osamstream&);

protected:
  virtual size_t xsgetn(isamstream&, char*, size_t);

private:
  char_buffer buffer;
  std::vector<char*> fields;
  bool reflist_open;
};

samio::samio() : buffer(32768), reflist_open(false) {
}

samio::samio(const char* text, std::streamsize textsize)
  : buffer(32768), reflist_open(false) {
  prepare_line_buffer(buffer, text, textsize);
}

samio::~samio() {
}

size_t samio::xsgetn(isamstream& stream, char* buffer, size_t length) {
  return rdbuf_sgetn(stream, buffer, length);
}

bool samio::get(isamstream& stream, collection& headers) {
  headers.clear();
  headers.refseqs_in_headers = true;
  set_cindex(headers);

  while (peek(buffer, stream) == '@') {
    int nfields = getline(buffer, stream, fields);
    headers.push_back(string(fields[0], fields[nfields] - fields[0] - 1),
		      add_header | add_refseq | add_refname);
  }

  // Further refsequences can be added if there were no @SQ headers.
  reflist_open = headers.ref_empty();
  headers.refseqs_in_headers = ! reflist_open;

  return ! (headers.empty() && peek(buffer, stream) == EOF);
}

bool samio::get(isamstream& stream, alignment& aln) {
  int nfields = getline(buffer, stream, fields);
  if (nfields <= 0)
    return false;

  //aln.p->h.cindex = header_cindex;
  // FIXME when looking up RNAME, MRNM, if reflist_open then can add unknown
  aln.assign(nfields, fields, header_cindex);
  return true;
}

void samio::flush(osamstream& stream) {
  while (buffer.size() > 0)
    buffer.begin += stream.rdbuf()->sputn(buffer.begin, buffer.size());

  buffer.clear();
}

void samio::put(osamstream& stream, const collection& headers) {
  for (collection::const_iterator it = headers.begin();
       it != headers.end(); ++it) {
    // FIXME Don't cons up a string
    string text = it->str();
    if (text.length() + 1 > buffer.available()) {
      flush(stream);
      buffer.reserve(text.length() + 1);
    }

    memcpy(buffer.end, text.data(), text.length());
    buffer.end += text.length();
    *buffer.end++ = '\n';
  }
}

void samio::put(osamstream& stream, const alignment& aln) {
  // TODO  If sync() ever affects more than just bin(), uncomment this!
  // aln.sync();

  size_t approx_length = aln.sam_length() + 1;
  if (approx_length > buffer.available()) {
    flush(stream);
    buffer.reserve(approx_length);
  }

  buffer.end = format_sam(buffer.end, aln, stream);
  *buffer.end++ = '\n';
}


// Gzipped SAM files
// =================

class gzsamio : public samio {
public:
  gzsamio(const char* text, std::streamsize textsize);
  virtual ~gzsamio();

protected:
  virtual size_t xsgetn(isamstream&, char*, size_t);
};

gzsamio::gzsamio(const char*, std::streamsize) {
  throw std::logic_error("gzsamio not implemented");
}

gzsamio::~gzsamio() {
}

size_t gzsamio::xsgetn(isamstream&, char*, size_t) {
  // TODO Read from rdbuf() and decompress
  throw std::logic_error("gzsamio::xsgetn() not implemented");
}

// *** virtual new

// Construct a new concrete sambamio by reading the first few bytes
// from the stream to determine what type of file it is.
sambamio* sambamio::new_in(isamstream& stream) {
  char buffer[BGZF::hsize];

  std::streamsize n = 0;
  while (n < BGZF::hsize && ! stream.eof())
    n += rdbuf_sgetn(stream, buffer + n, sizeof buffer - n);

  if (BGZF::is_bgzf_header(buffer, n))
    return new bamio(buffer, n);
  else if (BGZF::is_gzip_header(buffer, n))
    return new gzsamio(buffer, n);
  else
    return new samio(buffer, n);
}

// Construct a new concrete sambamio according to MODE.
sambamio* sambamio::new_out(std::ios::openmode mode) {
  if (mode & std::ios::binary)
    return new bamio(mode & compressed);
  else if (mode & compressed)
    return new gzsamio(NULL, 0);
  else
    return new samio();
}

} // namespace sam
