#include "lib/sambamio.h"

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
  bamio(const char* text, std::streamsize textsize);
  virtual ~bamio() { }

  virtual bool get_headers(isamstream&);
  virtual bool get(isamstream& strm, alignment& aln);
  virtual void put(osamstream& strm, const alignment& aln);
};

samstream_base::bamio::bamio(const char* text, std::streamsize textsize) {
}

bool samstream_base::bamio::get_headers(isamstream& strm) {
}

bool samstream_base::bamio::get(isamstream& strm, alignment& aln) {
}

void samstream_base::bamio::put(osamstream& strm, const alignment& aln) {
}


// *** samio

class samstream_base::samio : public samstream_base::sambamio {
public:
  samio(const char* text, std::streamsize textsize);
  virtual ~samio();

  virtual bool get_headers(isamstream&);
  virtual bool get(isamstream&, alignment&);
  virtual void put(osamstream&, const alignment&);

private:
  char* buffer;
  char* begin;
  char* end;
  size_t bufsize;

  std::vector<char*> fields;
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

bool samstream_base::samio::get_headers(isamstream& strm) {
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

void samstream_base::samio::put(osamstream& strm, const alignment& aln) {
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
samstream_base::sambamio::new_out(std::streambuf* sbuf, openmode mode) {
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
