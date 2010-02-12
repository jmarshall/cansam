#include "lib/sambamio.h"

#include <streambuf>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstddef>
#include <cstring>

#include <iostream> // FIXME NUKE-ME

#include <zlib.h>

#include "sam/alignment.h"
#include "sam/exception.h"
#include "sam/stream.h"
#include "lib/utilities.h"
#include "lib/wire.h"

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
enum { hsize = 18 };

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
inline int block_size(const char* text) {
  return convert::uint16(&text[16]) + 1;
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

  // Produce more available space beyond  end  by flushing, or, if the unread
  // characters are already at the start of the buffer, by enlarging it.
  // Updates begin/end and the specified pointers accordingly.
  char* flush_make_space(char* ptr, std::vector<char*>& ptrvec);

private:
  char* array;
  size_t capacity;
};

// Produces more space at the end of the buffer by moving the remaining
// unread characters in [begin,end) -- not including the sentinel -- to the
// start of the buffer, or, if  begin  is already at the start of the buffer,
// by enlarging it; rewires pointers, including  ptr  and those in  ptrvec,
// to point to the same places in the new location.  FIXME NUKE THIS COMMENT
char* char_buffer::flush_make_space(char* ptr, std::vector<char*>& ptrvec) {
  char* oldarray = NULL;
  char* oldbegin = begin;

  if (begin > array) {
    memmove(array, begin, end - begin);
  }
  else {
    size_t newcapacity = capacity * 2;
    char* newarray = new char[newcapacity];
    oldarray = array;
    memcpy(newarray, begin, end - begin);
    array = newarray;
    capacity = newcapacity;
  }

  begin = array;
  end = &begin[end - oldbegin];
  ptr = &begin[ptr - oldbegin];

  for (std::vector<char*>::iterator it = ptrvec.begin();
       it != ptrvec.end(); ++it)
    *it = &begin[*it - oldbegin];

  delete [] oldarray;

  return ptr;
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
	  if (b.available() == 0)  s = b.flush_make_space(s, fields);
	  *b.end = '\n';

	  *s++ = '\0';
	  fields.push_back(s);
	}

	break;
      }
      else {
	// This is the sentinel, and there are more characters to be read.
	s = b.flush_make_space(s, fields);

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
  virtual ~bamio();

  virtual bool get(isamstream&, collection&);
  virtual bool get(isamstream&, alignment&);
  virtual void put(osamstream&, const alignment&);

  virtual size_t xsgetn(isamstream&, char*, size_t);

private:
  size_t read(isamstream&, void*, size_t);
  int32_t read_int32(isamstream&);
  void read_refinfo(isamstream& stream, string& name, coord_t& length);

  bool underflow(isamstream&);
  void fill_cdata(isamstream&);
  size_t inflate_into_buffer(char*, size_t);

  char_buffer buffer;
  char_buffer cdata;

  size_t header_text_length;
};

bamio::bamio(const char* text, std::streamsize textsize)
  : buffer(65536), cdata(65536) {
  memcpy(cdata.end, text, textsize);
  cdata.end += textsize;
}

bamio::~bamio() {
}

// Fill  cdata  by reading from the streambuf.  Reads as much as will fit into
// the buffer, which is probably several BGZF blocks.  This reduces the number
// of read(2) system calls for sequential access, and (with modern disk block
// sizes) shouldn't mean too much data is discarded by any subsequent seek().
void bamio::fill_cdata(isamstream& stream) {
  cdata.flush();
  cdata.end += rdbuf_sgetn(stream, cdata.end, cdata.available());
}

// FIXME Maybe don't have underflow() skip the header so we can do
// inflateInit() normal not raw mode and get CRC checking

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
size_t bamio::inflate_into_buffer(char* data, size_t length) {
  z_stream z;
  z.zalloc = Z_NULL;
  z.zfree  = Z_NULL;

  z.next_in  = reinterpret_cast<unsigned char*>(data);
  z.avail_in = length;
  if (inflateInit2(&z, -15) != Z_OK)
    throw bad_format(zlib_error_message("inflateInit2", z));

  buffer.clear();
  z.next_out  = reinterpret_cast<unsigned char*>(buffer.end);
  z.avail_out = buffer.available();
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

// Decompress one BGZF block into  buffer,  which is assumed to be previously
// empty, from  cdata,  which is refilled by reading from the streambuf if
// necessary.  Returns true if  buffer  is nonempty afterwards.
bool bamio::underflow(isamstream& stream) {
  if (cdata.size() < BGZF::hsize) {
    fill_cdata(stream);
    // If there's still no data, we're cleanly at EOF.
    if (cdata.size() == 0)  return false;
  }

  if (BGZF::is_bgzf_header(cdata.begin, cdata.size())) {
    size_t blockonly_length = BGZF::block_size(cdata.begin) - BGZF::hsize;
    cdata.begin += BGZF::hsize;
    if (cdata.size() < blockonly_length) {
      fill_cdata(stream);
      if (cdata.size() < blockonly_length)
	throw bad_format(make_string()
	    << "Truncated BGZF block (expected " << blockonly_length
	    << " bytes after header; got " << cdata.size() << ")");
    }

#if 1
    cdata.begin += inflate_into_buffer(cdata.begin, blockonly_length);
    cdata.begin += 8; // FIXME skip the GZIP member footer... should check it
#else
    size_t n;
    cdata.begin += n = inflate_into_buffer(cdata.begin, blockonly_length);
    //std::clog << "inflated: blockonly_length:" << blockonly_length << "; returned:" << n << "\n";
    cdata.begin += 8; // FIXME skip the GZIP member footer... should check it
#endif
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
  set_cindex(headers);

  header_text_length = read_int32(stream);

  char_buffer header_text_buffer(32768);
  prepare_line_buffer(header_text_buffer);
  std::vector<char*> fields;
  while (peek(header_text_buffer, stream) == '@') {
    int nfields = getline(header_text_buffer, stream, fields);
    headers.push_back(string(fields[0], fields[nfields] - fields[0] - 1),
		      add_header | add_refname);
  }

  if (headers.refnames.empty()) {
#if 0
    string name;
    coord_t length;
    int ref_count = read_int32(stream);
#endif

    std::clog << "no @SQ headers\n";
  }
  else {
    string name;
    coord_t length;
    int ref_count = read_int32(stream);
    for (int index = 0; index < ref_count; index++) {
      read_refinfo(stream, name, length);
      // FIXME more checking...
      refsequence* rhdr = headers.refnames[name];
      //rhdr->index_ = index;
      headers.refseqs.push_back(rhdr);
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

void bamio::put(osamstream& stream, const alignment& aln) {
  aln.sync();

  size_t length = min(aln.p->size(), buffer.available());
  memcpy(buffer.end, aln.p->data(), length);
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

  // FIXME emit the rest

#if 0
  commit(stream);
#endif
  #warning bamio::put(aln) unimplemented
  stream.eof();
}


// Text SAM files
// ==============

class samio : public sambamio {
public:
  samio(const char* text, std::streamsize textsize);
  virtual ~samio();

  virtual bool get(isamstream&, collection&);
  virtual bool get(isamstream&, alignment&);
  virtual void put(osamstream&, const alignment&);

  virtual size_t xsgetn(isamstream&, char*, size_t);

protected:
  samio() : buffer(32768), reflist_open(false) { }

private:
  char_buffer buffer;
  std::vector<char*> fields;
  bool reflist_open;
};

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
  set_cindex(headers);

  while (peek(buffer, stream) == '@') {
    int nfields = getline(buffer, stream, fields);
    headers.push_back(string(fields[0], fields[nfields] - fields[0] - 1),
		      add_header | add_refseq | add_refname);
  }

  // Further refsequences can be added if there were no @SQ headers.
  reflist_open = headers.refseqs.empty();

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

void samio::put(osamstream& stream, const alignment& aln) {
  aln.sync();

  #warning samio::put(aln) unimplemented
  stream.eof();
}


// Gzipped SAM files
// =================

class gzsamio : public samio {
public:
  gzsamio(const char* text, std::streamsize textsize);
  virtual ~gzsamio();

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
sambamio* sambamio::new_in(std::streambuf* sbuf) {
  char buffer[BGZF::hsize];
  std::streamsize n = sbuf->sgetn(buffer, sizeof buffer);
  // FIXME eofbit?

  if (BGZF::is_bgzf_header(buffer, n))
    return new bamio(buffer, n);
  else if (BGZF::is_gzip_header(buffer, n))
    return new gzsamio(buffer, n);
  else
    return new samio(buffer, n);
}

// Construct a new concrete sambamio according to MODE.
sambamio* sambamio::new_out(std::streambuf* /*sbuf*/, std::ios::openmode mode) {
  if (mode & std::ios::binary)
    return new bamio(NULL, 0 /*, mode & compressed*/);
  else if (mode & compressed)
    return new gzsamio(NULL, 0);
  else
    return new samio(NULL, 0);
}

} // namespace sam
