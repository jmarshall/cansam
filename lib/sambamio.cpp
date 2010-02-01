#include "lib/sambamio.h"

#include <streambuf>
#include <string>
#include <vector>

#include <cstring>

#include <zlib.h>

#include "sam/alignment.h"
#include "sam/collection.h"
#include "sam/exception.h"
#include "sam/stream.h"
#include "lib/bgzf.h"
#include "lib/buffer.h"
#include "lib/utilities.h"
#include "lib/wire.h"

using std::string;

namespace sam {

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

// *** bamio

class samstream_base::bamio : public samstream_base::sambamio {
public:
  bamio(const char* text, std::streamsize textsize);
  virtual ~bamio() { }

  virtual bool get(isamstream&, collection&);
  virtual bool get(isamstream&, alignment&);
  virtual void put(osamstream&, const alignment&);

private:
  size_t read(isamstream&, void*, size_t);
  bool underflow(isamstream&);
  void fill_cdata(isamstream&);
  size_t inflate_into_buffer(char*, size_t);

  read_buffer buffer;
  read_buffer cdata;
};

samstream_base::bamio::bamio(const char* text, std::streamsize textsize)
  : buffer(65536), cdata(65536) {
  memcpy(cdata.end, text, textsize);
  cdata.end += textsize;
}

// Fill  cdata  by reading from the streambuf.  Reads as much as will fit into
// the buffer, which is probably several BGZF blocks.  This reduces the number
// of read(2) system calls for sequential access, and (with modern disk block
// sizes) shouldn't mean too much data is discarded by any subsequent seek().
void samstream_base::bamio::fill_cdata(isamstream& stream) {
  cdata.flush();
  // FIXME Probably remove rdbuf_sgetn and put its code here instead.
  cdata.end += stream.rdbuf_sgetn(cdata.end, cdata.capacity() - cdata.size());
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
size_t samstream_base::bamio::inflate_into_buffer(char* data, size_t length) {
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

// Decompress one BGZF block into  buffer,  which is assumed to be previously
// empty, from  cdata,  which is refilled by reading from the streambuf if
// necessary.  Returns true if  buffer  is nonempty afterwards.
bool samstream_base::bamio::underflow(isamstream& stream) {
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

    cdata.begin += inflate_into_buffer(cdata.begin, blockonly_length);
    return true;
  }
  else
    throw bad_format("Invalid BGZF block header");
}

size_t samstream_base::bamio::read(isamstream& stream,
				   void* destv, size_t desired_length) {
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

bool samstream_base::bamio::get(isamstream& stream, collection& headers) {
  /* read magic
     read l_text, text
     headers.clear();
     for each header line, headers.push_back(line);
     nref = read n_ref
     for (
  */
  #warning bamio::get(headers) unimplemented
  stream.eof();
  headers.begin();
  return false;
}

bool samstream_base::bamio::get(isamstream& stream, alignment& aln) {
  uint32_t rest_length;

  size_t n = read(stream, &rest_length, sizeof rest_length);
  if (n < sizeof rest_length) {
    if (n == 0)  return false;
    else  throw bad_format("Truncated BAM alignment record");
  }

  convert::set_uint32(rest_length);
  aln.reserve(rest_length); // FIXME er plus/minus a bit

  n = read(stream, &aln.p->c.rindex, rest_length);
  if (n < rest_length)
    throw bad_format(make_string()
	<< "Truncated BAM alignment record (got " << n
	<< " bytes of an expected remainder of " << rest_length << ")");

  aln.p->h.cindex = 5; // FIXME stream.cindex;
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

void samstream_base::bamio::put(osamstream& stream, const alignment& aln) {
  #warning bamio::put(aln) unimplemented
  stream.eof();
  aln.find("X0");
}


// *** samio

class samstream_base::samio : public samstream_base::sambamio {
public:
  samio(const char* text, std::streamsize textsize);
  virtual ~samio();

  virtual bool get(isamstream&, collection&);
  virtual bool get(isamstream&, alignment&);
  virtual void put(osamstream&, const alignment&);

private:
  char* buffer;
  char* begin;
  char* end;
  size_t bufsize;

  std::vector<char*> fields;
  int peek(isamstream&);
  int getline(isamstream&);
  char* flush_buffer(char*);
};

samstream_base::samio::samio(const char* text, std::streamsize textsize) {
  bufsize = 32768;
  buffer = new char[bufsize];

  end = &buffer[bufsize - 1];
  begin = end - textsize;
  memcpy(begin, text, textsize);
  *end = '\n';
}

samstream_base::samio::~samio() {
  delete [] buffer;
}

/* Pointers to the read buffer are arranged as follows:

  [---------ABCDEF/GHIJKL/MNO
  ^         ^                ^
  buffer    begin            end */

// Produces more space at the end of the buffer by moving the remaining
// unread characters in [begin,end) -- not including the sentinel -- to the
// start of the buffer, or, if  begin  is already at the start of the buffer,
// by enlarging it; rewires pointers, including  ptr  and those in  fields,
// to point to the same places in the new location.
char* samstream_base::samio::flush_buffer(char* ptr) {
  char* oldbuffer = NULL;
  char* oldbegin = begin;

  if (begin > buffer) {
    memmove(buffer, begin, end - begin);
  }
  else {
    bufsize *= 2;
    char* newbuffer = new char[bufsize];
    oldbuffer = buffer;
    memcpy(newbuffer, begin, end - begin);
    buffer = newbuffer;
  }

  begin = buffer;
  end = &begin[end - oldbegin];
  ptr = &begin[ptr - oldbegin];

  for (std::vector<char*>::iterator it = fields.begin();
       it != fields.end(); ++it)
    *it = &begin[*it - oldbegin];

  delete [] oldbuffer;

  return ptr;
}

int samstream_base::samio::peek(isamstream& stream) {
  #warning samio::peek() unimplemented
  #warning - or is ungetline() a better idea?
  stream.eof();
  return EOF;
}

/* Reads a newline-terminated line of tab-delimited text into  fields,
and returns the number of fields present (or 0 at EOF).

*/
int samstream_base::samio::getline(isamstream& stream) {
  fields.clear();
  fields.push_back(begin);

  char* s = begin;
  while (true)
    if (*s == '\t') {
      *s++ = '\0';
      fields.push_back(s);
    }
    else if (*s == '\n') {
      if (s < end) {
	// This is a real \n, so a properly-terminated line has been read.
	if (s > begin && s[-1] == '\r')  s[-1] = '\0', fields.push_back(s++);
	else  *s++ = '\0', fields.push_back(s);

	break;
      }
      else if (stream.eof()) {
	// This is the sentinel and no further characters are forthcoming.
	// If any characters have been read, they constitute a final line
	// (which is unterminated); otherwise we are properly at EOF.
	if (s > begin) {
	  // Move the sentinel one character later, to make room for a \0.
	  end++;
	  if (end >= &buffer[bufsize])  s = flush_buffer(s);
	  *end = '\n';

	  *s++ = '\0';
	  fields.push_back(s);
	}

	break;
      }
      else {
	// This is the sentinel, and there are more characters to be read.
	s = flush_buffer(s);

	// Buffer space available to be filled is [end,bufsize-1), leaving
	// one character to spare for the sentinel.
	std::streamsize n = stream.rdbuf()->sgetn(end, &buffer[bufsize-1]-end);
	if (n == 0) {
	  if (stream.setstate_wouldthrow(eofbit))
	    throw eof_exception();
	}

	end += n;
	*end = '\n';
      }
    }
    else
      s++;

  begin = s;
  return fields.size() - 1;
}

bool samstream_base::samio::get(isamstream& stream, collection& headers) {
  headers.clear();

  while (peek(stream) == '@') {
    int nfields = getline(stream);
    headers.push_back(string(fields[0], fields[nfields] - fields[0]));
  }

  #warning samio::get(headers) unimplemented
  return false;
}

bool samstream_base::samio::get(isamstream& stream, alignment& aln) {
  int nfields = getline(stream);
#if 0
std::clog << "getline returned " << nfields << ":";
for (int i = 0; i < nfields; i++)
  std::clog << " {" << fields[i] << "}";
std::clog << "\n";
#endif
  if (nfields <= 0)
    return false;

  aln.assign(nfields, fields /*, stream.refthingie*/);
  return true;
}

void samstream_base::samio::put(osamstream& stream, const alignment& aln) {
  #warning samio::put(aln) unimplemented
  stream.eof();
  aln.find("X0");
}


// *** virtual new

// Construct a new concrete sambamio by reading the first few bytes
// from the stream to determine what type of file it is.
samstream_base::sambamio*
samstream_base::sambamio::new_in(std::streambuf* sbuf) {
  // BAM files begin with an RFC1952 GZIP member header with a BAM-specific
  // "BC" extra subfield, typically as follows:
  //
  //  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17   Offset from start
  // 1f 8b .. *4 .. .. .. .. .. .. 06 00 42 43 02 00 .. ..   BAM signature
  // magic    ^bit 2 is FEXTRA flag      ^B ^C

  char buf[16];
  std::streamsize n = sbuf->sgetn(buf, sizeof buf);
  // FIXME eofbit?

  // Does the stream start with the GZIP magic number?
  if (n >= 2 && buf[0] == '\x1f' && buf[1] == '\x8b') {
    // If so, does it also have a "BC" extra subfield?
    // FIXME Allow for other extra subfields in addition to the "BC" one.
    // FIXME Theoretically should allow for other extra subfields.
    if (n >= 16 && (buf[3] & 4) && memcmp(&buf[10], "\6\0\x42\x43\2\0", 6) == 0)
      return new bamio(buf, n);
    else
#if 0
      return new gzsamio(buf, n);
#else
      throw ".sam.gz not yet implemented";
#endif
  }
  else
    return new samio(buf, n);
}

// Construct a new concrete sambamio according to MODE.
samstream_base::sambamio*
samstream_base::sambamio::new_out(std::streambuf* /*sbuf*/, openmode mode) {
  if (mode & binary)
    return new bamio(NULL, 0 /*, mode & compressed*/);
#if 0
  else if (mode & compressed)
    return new gzsamio(NULL, 0);
#endif
  else
    return new samio(NULL, 0);
}

} // namespace sam
