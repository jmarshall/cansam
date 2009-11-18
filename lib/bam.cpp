#include <streambuf>
#include <vector>

#include <cstring>

#include "sam/stream.h"
#include "lib/zio.h"
#include "lib/wire.h"

namespace sam {

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
  bamio(const char* buffer, std::streamsize bufsize);
  virtual ~bamio() { }

  virtual bool get_headers(samstream_base&);
  virtual bool get(samstream_base& strm, alignment& aln);
};

bool samstream_base::bamio::get_headers(samstream_base& strm) {
}

bool samstream_base::bamio::get(samstream_base& strm, alignment& aln) {
}

// *** samio

class samstream_base::samio : public samstream_base::sambamio {
public:
  samio(const char* buffer, std::streamsize bufsize);
  virtual ~samio();

  virtual bool get_headers(samstream_base&);
  virtual bool get(samstream_base& strm, alignment& aln);

private:
  char buffer[65537];
  char* begin;
  char* end;

  std::vector<char*> fields;
  bool getline();
};

samstream_base::samio::samio(const char* buffer, std::streamsize bufsize) {
}

samstream_base::samio::~samio() {
}

bool samstream_base::samio::getline() {
  if (begin >= end)
    return false;

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
	if (s[-1] == '\r') {
	  // ...
	}
	break;
      }
      else {
	// move, realloc, redo pointers
      }
    }
    else
      s++;
 
}

bool samstream_base::samio::get(samstream_base& stream, alignment& aln) {
  // An alignment in SAM format is a tab-separated line containing fields
  // ordered as:
  // qname flag rname pos mapq cigar mrname mpos isize seq qual aux...
  // 0     1    2     3   4    5     6      7    8     9   10   11...


  //gettabbedline(fields);
  // aln.assign(fields, something);
}

// *** gzsamio

class samstream_base::gzsamio : public samstream_base::samio {
public:
  gzsamio(const char* buffer, std::streamsize bufsize);
  virtual ~gzsamio() { }

  virtual bool get_headers(samstream_base&);
  virtual bool get(samstream_base& strm, alignment& aln);
};

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

  // Does the stream start with the GZIP magic number?
  if (n >= 2 && buf[0] == '\x1f' && buf[1] == '\x8b') {
    // If so, does it also have a "BC" extra subfield?
    // FIXME Allow for other extra subfields in addition to the "BC" one.
    // FIXME Theoretically should allow for other extra subfields.
    if (n >= 16 && (buf[3] & 4) && memcmp(&buf[10], "\6\0\x42\x43\2\0", 6) == 0)
      return new bamio(buf, n);
    else
      return new gzsamio(buf, n);
  }
  else
    return new samio(buf, n);
}

// Construct a new concrete sambamio according to MODE.
samstream_base::sambamio*
samstream_base::sambamio::new_out(std::streambuf*, openmode mode) {
  if (mode & binary)
    return new bamio(NULL, 0 /*, mode & compressed*/);
  else if (mode & compressed)
    return new gzsamio(NULL, 0);
  else
    return new samio(NULL, 0);
}

} // namespace sam
