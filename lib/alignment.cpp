#define __STDC_LIMIT_MACROS
#include <stdint.h>
#if !(defined UINT8_MAX && defined UINT16_MAX && defined INT32_MAX)
#error - This implementation requires 8-, 16-, and 32-bit integer types
#endif

#include "sam/alignment.h"

#include <stdexcept>
#include <string>
#include <cstring>
#include <cstdlib>  // for strtol()
#include <climits>  // for UCHAR_MAX

#include <iostream> // FIXME NUKE-ME

#include "sam/collection.h"
#include "sam/exception.h"

#include "lib/utilities.h"
#include "lib/wire.h"

using std::string;

namespace sam {

// 1. Accessors
//=============

// For order(); 0, FIRST_IN_PAIR, SECOND_IN_PAIR, FIRST_IN_PAIR|SECOND_IN_PAIR.
const int alignment::order_value[4] = { 0, -1, +1, 0 };

string alignment::cigar() const {
  return "FIXME"; // FIXME
}

scoord_t alignment::cigar_span() const {
  return 0;  // FIXME
}

void alignment::pack_seq(char* dest, const char* seq, int seq_length) {
  static char encode[UCHAR_MAX + 1];
  typedef unsigned char uchar;

  if (encode[uchar('A')] == '\0') {
    for (int i = 0; i <= UCHAR_MAX; i++)  encode[i] = 16;
    encode[uchar('=')] = 0;
    encode[uchar('A')] = encode[uchar('a')] = 1;
    encode[uchar('C')] = encode[uchar('c')] = 2;
    encode[uchar('M')] = encode[uchar('m')] = 3;
    encode[uchar('G')] = encode[uchar('g')] = 4;
    encode[uchar('R')] = encode[uchar('r')] = 5;
    encode[uchar('S')] = encode[uchar('s')] = 6;
    encode[uchar('V')] = encode[uchar('v')] = 7;
    encode[uchar('T')] = encode[uchar('t')] = 8;
    encode[uchar('W')] = encode[uchar('w')] = 9;
    encode[uchar('Y')] = encode[uchar('y')] = 10;
    encode[uchar('H')] = encode[uchar('h')] = 11;
    encode[uchar('K')] = encode[uchar('k')] = 12;
    encode[uchar('D')] = encode[uchar('d')] = 13;
    encode[uchar('B')] = encode[uchar('b')] = 14;
    encode[uchar('N')] = encode[uchar('n')] = 15;
    encode[uchar('.')] = 15;
  }

  const unsigned char* unpacked = reinterpret_cast<const unsigned char*>(seq);

  int even_length = seq_length & ~1;
  for (int i = 0; i < even_length; i += 2) {
    char hi = encode[*unpacked++];
    char lo = encode[*unpacked++];
    if ((hi | lo) & 16)
      throw bad_format(make_string()
	  << "Invalid character ('" << unpacked[(hi & 16)? -2 : -1]
	  << "') in sequence string");
    *dest++ = (hi << 4) | lo;
  }

  if (even_length < seq_length) {
    char hi = encode[*unpacked];
    if (hi & 16)
      throw bad_format(make_string()
	  << "Invalid character ('" << *unpacked << "') in sequence string");
    *dest = (hi << 4);
  }
}

// Unpack the BAM-encoded sequence data, which is of course not necessarily
// NUL-terminated.  Assumes that DEST is already of sufficient length for
// all this (c.f. unpack_seq()).
void alignment::unpack_seq(string::iterator dest,
			   const char* raw_seq, int seq_length) {
  static const char decode[] =
    "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
    "A=AAACAMAGARASAVATAWAYAHAKADABAN"  //  1 = A
    "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"  //  2 = C
    "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
    "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"  //  4 = G
    "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
    "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
    "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
    "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"  //  8 = T
    "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
    "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"  // ...and similarly for columns w.r.t.
    "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"  // the second character of each pair.
    "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
    "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
    "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
    "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN"; // 15 = N

  const unsigned char* packed = reinterpret_cast<const unsigned char*>(raw_seq);

  int even_length = seq_length & ~1;
  for (int i = 0; i < even_length; i += 2) {
    int ndx = 2 * *packed++;
    *dest++ = decode[ndx];
    *dest++ = decode[ndx+1];
  }

  if (even_length < seq_length) {
    int ndx = 2 * *packed;
    *dest = decode[ndx];
  }
}

void alignment::pack_qual(char* dest, const char* qual, int seq_length) {
  for (int i = 0; i < seq_length; i++) {
    char q = *qual++;
    if (q < 33 || q > 126)
      throw bad_format(make_string()
	  << "Invalid character ('" << q << "') in quality string");
    *dest++ = q - 33;
  }
}

void alignment::unpack_qual(char* dest, const char* phred, int seq_length) {
  for (int i = 0; i < seq_length; i++) {
    char q = *phred++;
    if (q < 0)  q = 0;
    else if (q > 126 - 33)  q = 126 - 33;

    *dest++ = q + 33;
  }
}

bool operator< (const alignment& a, const alignment& b) {
  // Treating these as unsigned means -1 (unmapped) sorts last.
  unsigned a_rindex = unsigned(a.rindex());
  unsigned b_rindex = unsigned(b.rindex());
  if (a_rindex != b_rindex)  return (a_rindex < b_rindex);

  if (a.pos() != b.pos())  return a.pos() < b.pos();

  int dqname = strcmp(a.qname_c_str(), b.qname_c_str());
  if (dqname != 0)  return (dqname < 0);

  return a.order() < b.order();
}

// 1a. Mutators
//=============

void alignment::set_qname(const char* qname_data, int qname_length) {
  char* qbuffer =
    replace_gap(p->data() + p->name_offset(), p->data() + p->cigar_offset(),
		qname_length + 1);

  p->c.name_length = qname_length + 1;
  memcpy(qbuffer, qname_data, qname_length);
  qbuffer[qname_length] = '\0';
}

void alignment::set_flags(int flags) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.flags = flags;
}

// void alignment::set_rindex(int rindex)
// void alignment::set_rname(const std::string& rname)

void alignment::set_pos(coord_t pos) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.zpos = pos - 1;
  p->c.bin = unknown_bin;
}

void alignment::set_zpos(coord_t zpos) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.zpos = zpos;
  p->c.bin = unknown_bin;
}

void alignment::set_mapq(int mapq) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.mapq = mapq;
}

// void alignment::set_cigar(const std::string& cigar)
// void alignment::set_mate_rindex(int mate_rindex)
// void alignment::set_mate_rname(const std::string& mate_rname)

void alignment::set_mate_pos(coord_t pos) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.mate_zpos = pos - 1;
}

void alignment::set_mate_zpos(coord_t zpos) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.mate_zpos = zpos;
}

void alignment::set_isize(scoord_t isize) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.isize = isize;
}

// 2. Other biostuffs
//===================

int calc_zbin(coord_t zpos, coord_t zright) {
  if (zpos >> 14 == zright >> 14)  return ((1 << 15) - 1) / 7 + (zpos >> 14);
  if (zpos >> 17 == zright >> 17)  return ((1 << 12) - 1) / 7 + (zpos >> 17);
  if (zpos >> 20 == zright >> 20)  return ((1 <<  9) - 1) / 7 + (zpos >> 20);
  if (zpos >> 23 == zright >> 23)  return ((1 <<  6) - 1) / 7 + (zpos >> 23);
  if (zpos >> 26 == zright >> 26)  return ((1 <<  3) - 1) / 7 + (zpos >> 26);
  return 0;
}

int calc_bin(coord_t pos, coord_t right) {
  return calc_zbin(pos - 1, right - 1);
}

void alignment::sync() const {
  bin();
}

// 3. Infrastructure
//==================

/*+---------------+---------+--...--+-...-+-...--+--...--+------...------+
  | capacity, etc | bamcore | cigar | seq | qual | rname | aux fields... |
  +---------------+---------+--...--+-...-+-...--+--...--+------...------+*/

/* An alignment::block representing an uninitialised alignment, which will be
shared by all the default-constructed alignments.  (It's shared so that we can
have a constant-time default constructor.)  This block lies about its capacity
so that tests of the form "p->capacity() < some_size" always trigger when  p
is the empty block.  */
struct alignment::block alignment::empty_block = {
  { 0 /* 41, if truth be told */, 0 },
  { 33, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0 },
  { '\0' /* an empty qname C-string */ }
};

// Allocate a new block and initialise its capacity field.
alignment::block* alignment::block::create(int capacity) {
  char* cp = new char[capacity];
  block* p = reinterpret_cast<block*>(cp);

  p->h.capacity = capacity;
  return p;
}

// Copy the block contents (but don't overwrite DEST's capacity field),
// assuming that the destination's capacity suffices for the source's size.
void alignment::block::copy(block* dest, const block* src) {
  memcpy(&dest->h.cindex, &src->h.cindex,
	 src->size() - sizeof(src->h.capacity));
}

// Deallocate the block.
void alignment::block::destroy(block* p) {
  // FIXME do if != &empty_block test here?  Otherwise inline?
  char* cp = reinterpret_cast<char*>(p);
  delete [] cp;
}

/* Resize the alignment's block, unsharing it if it is currently the shared
empty block, while maintaining the existing contents.  Used by mutators to
ensure the block is sufficiently sized and writable.  */
void alignment::resize_unshare_copy(int size) {
  block* newp = block::create(size);
  block* oldp = p;
  block::copy(newp, oldp);
  p = newp;
  if (oldp != &empty_block)  block::destroy(oldp);
}

/* Resize the alignment's block, unsharing it if it is currently the shared
empty block, but not maintaining the existing contents.  Used by assignments
to ensure the block is sufficiently sized and writable.  */
void alignment::resize_unshare_discard(int size) {
  block* newp = block::create(size);
  block* oldp = p;
  p = newp;
  if (oldp != &empty_block)  block::destroy(oldp);
}

/* This assignment operator simply copies the pointed-to block, but must first
ensure that the destination is large enough and is not the shared empty block,
and must explicitly check for the self-assignment case.  */
alignment& alignment::operator= (const alignment& aln) {
  int newsize = aln.p->size();
  if (p->capacity() < newsize)
    resize_unshare_discard(newsize);
  else if (p == aln.p)
    return *this;

  block::copy(p, aln.p);
  return *this;
}

alignment::alignment(const alignment& aln)
  : p(block::create(aln.p->size())) {
  block::copy(p, aln.p);
}

#if 1
alignment& alignment::assign(const string& /*line*/) {
  // FIXME
  return *this;
}
#endif

// FIXME Where do these functions go?
int to_flags(const char* s) {
  int val;

  if (isdigit(s[0])) {
    char* slim;
    val = strtol(s, &slim, 0);
    if (*slim != '\0')
      throw bad_format("Trailing gunge in flags field");
  }
  else {
    // FIXME One day there'll be a defined text representation
    val = 0;
  }

  return val;
}

// FIXME  As in, the number of operators
int calc_cigar_length(const char* s) {
  int n = 0;

  if (s[0] == '*' && s[1] == '\0') {
    // An unspecified cigar string, written as "*", is represented as
    // empty in BAM, so leave  n  as 0.
  }
  else {
    while (*s)
      if (!isdigit(*s++))  n++;
  }

  return n;
}

#if 0
int sam_aux_length(char type, const char* value, int value_length) {
  switch (type) {
  case 'A':
    if (value_length != 1)
      throw bad_format("Type 'A' aux field has length other than 1");
    return 1;

  case 'i':
  case 'f':
  case 'd':
  case 'Z':
  case 'H':
    break;
  }
}
#else
// FIXME
int sam_aux_length(char, const char*, int) {
  return 20;
}
#endif

int lookup_in_refthingie(const char*) { return -1; }
int itoa(const char*) { return 0; }

int alignment::approx_sam_record_length() const {
  int len = 0;
  len += p->c.name_length;  // This includes the \0, so accounts for the tab
  len += 4 + 1;
  len += 5; // FIXME refname
  len += 10 + 1;
  len += 3 + 1;
  len += 5; // FIXME cigar

  return len;
}

template <typename IntType>
char* format_hex(char* dest, IntType val) {
  *dest++ = '0';
  *dest++ = 'x';

  IntType n = val;
  do { dest++; n >>= 4; } while (n != 0);

  char* destlim = dest;
  do { *--dest = "0123456789ABCDEF"[val & 0xf]; val >>= 4; } while (val != 0);

  return destlim;
}

void alignment::sam_record(char* dest, int /*dest_length*/) const {
#if 0
  if (! length_is_okay)
    return excess_needed;
#endif

  strcpy(dest, qname_c_str()), dest += p->c.name_length - 1;

  *dest++ = '\t';
  if (0 /*fmtflags & boolalpha*/)
    strcpy(dest, "FIXME"), dest += 5; // FIXME
  else if (0 /*fmtflags & hex*/)
    dest = format_hex(dest, flags());
  else
    dest = format::decimal(dest, flags());

  *dest++ = '\t';
  if (rindex() < 0)  *dest++ = '*';
  else  strcpy(dest, "REFNAME"), dest += 7; // FIXME

  *dest++ = '\t';
  dest = format::decimal(dest, pos());

  *dest++ = '\t';
  dest = format::decimal(dest, mapq());

  *dest++ = '\t';
  strcpy(dest, "CIGAR"), dest += 5; // FIXME

  *dest++ = '\t';
  if (mate_rindex() < 0)  *dest++ = '*';
  else if (mate_rindex() == rindex())  *dest++ = '=';
  else  strcpy(dest, "MATEREFNAME"), dest += 11; // FIXME

  *dest++ = '\t';
  // FIXME Do we need to do anything special for 0 or -1 i.e. unmapped?
  dest = format::decimal(dest, mate_pos());

  *dest++ = '\t';
  dest = format::decimal(dest, isize());

  for (const_iterator it = begin(); it != end(); ++it) {
    *dest++ = '\t';
    *dest++ = it->tag_[0], *dest++ = it->tag_[1];
    *dest++ = ':';
    *dest++ = it->type_; // FIXME sanitize BAM-only types
    *dest++ = ':';
    // FIXME Append the aux value
    strcpy(dest, "FIXME"), dest += 5;
  }
}

void alignment::assign(int nfields, const std::vector<char*>& fields /*, refthingie*/) {
  // An alignment in SAM format is a tab-separated line containing fields
  // ordered as:
  // qname flag rname pos mapq cigar mrname mpos isize seq qual aux...
  // 0     1    2     3   4    5     6      7    8     9   10   11...
  enum { qname, flag, rname, pos, mapq, cigar, mrname, mpos, isize, seq, qual,
	 firstaux };

  if (nfields < 11)
    throw bad_format("Too few fields in SAM record");

  int name_length = fields[qname+1] - fields[qname]; // including terminator
  int cigar_length = calc_cigar_length(fields[cigar]);
  int seq_length = fields[seq+1] - fields[seq] - 1;
  int qual_length = fields[qual+1] - fields[qual] - 1;

  if (seq_length == 1 && fields[seq][0] == '*')
    seq_length = 0;

  if (qual_length == 1 && fields[qual][0] == '*')
    qual_length = 0;
  else if (qual_length != seq_length) {
    if (seq_length == 0)
      throw bad_format("QUAL specified when SEQ is absent");
    else
      throw bad_format("SEQ and QUAL differ in length");
  }

  int size = sizeof(block_header) + sizeof(bamcore)
	     + name_length + (cigar_length * sizeof(uint32_t))
	     + (seq_length+1)/2 + seq_length;

  for (int i = firstaux; i < nfields; i++) {
    int length = fields[i+1] - fields[i] - 1;
    char* aux = fields[i];
    if (length < 5 || aux[2] != ':' || aux[4] != ':')
      throw bad_format("Malformatted aux field");
    size += 3 + sam_aux_length(aux[3], &aux[5], length - 5);
  }

  if (p->capacity() < size)
    resize_unshare_discard(size);

  p->h.cindex = 0; // FIXME somethingie

  p->c.rest_length = size - sizeof(block_header) - sizeof(p->c.rest_length);
  p->c.rindex = lookup_in_refthingie(fields[rname]); // a name or "*"
  p->c.zpos = itoa(fields[pos]) - 1; // a 1-based pos or 0
  p->c.name_length = name_length;
  p->c.mapq = itoa(fields[mapq]); // 0..255
  p->c.bin = -1;
  p->c.cigar_length = cigar_length;
  p->c.flags = to_flags(fields[flag]);
  p->c.read_length = seq_length;
  p->c.mate_rindex = lookup_in_refthingie(fields[mrname]); // a name or "*"
  p->c.mate_zpos = itoa(fields[mpos]) - 1; // a 1-based pos or 0
  p->c.isize = itoa(fields[isize]); // an insert size or 0

  // In BAM: int32_t read_len is given, seq is (read_len+1)/2 bytes;
  //         qual is read_len bytes (Phred qualities), maybe 0xFF x read_len.

  // In SAM: either seq is "*" so read_len is 0, or seq is read_len base chars;
  //         qual is "*" (so 0xFF x read_len) or it's P-33 x read_len qualchars.

  // cigar

  pack_seq(p->data() + p->seq_offset(), fields[seq], seq_length);
  if (qual_length > 0)
    pack_qual(p->data() + p->qual_offset(), fields[qual], seq_length);
  else
    memset(p->data() + p->qual_offset(), '\xff', seq_length);

  memcpy(p->data() + p->name_offset(), fields[qname], name_length);

  for (int i = firstaux; i < nfields; i++) {
    char* aux = fields[i];

    switch (aux[3]) {
    case 'A':
      throw std::logic_error("Aux 'A' field not implemented");  // TODO

    case 'i':
      push_back(aux, 37); //to_int(&aux[5]));
      break;

    case 'f':
      throw std::logic_error("Aux 'f' field not implemented");  // TODO

    case 'd':
      throw std::logic_error("Aux 'd' field not implemented");  // TODO

    case 'Z':
      push_back(aux, &aux[5]);
      break;

    case 'H':
      throw std::logic_error("Aux 'H' field not implemented");  // TODO

    default:
      throw bad_format(make_string()
	  << "Aux field '" << aux[0] << aux[1] << "' has invalid type ('"
	  << aux[3] << "') for SAM format");
    }
  }
}

#if 0
int alignment::to_flags(const string& str) {
  const char* s = str.data();
  return to_flags(s, s + str.length());
}
#endif

// 4. Iterators
//=============

int alignment::tagfield::size() const {
  switch (type_) {
  case 'A':
    return 2 + 1 + 1;

  case 'c':
  case 'C':
    return 2 + 1 + sizeof(int8_t);

  case 's':
  case 'S':
    return 2 + 1 + sizeof(int16_t);

  case 'i':
  case 'I':
    return 2 + 1 + sizeof(int32_t);

  case 'f':
    return 2 + 1 + 4;

  case 'd':
    return 2 + 1 + 8;

  case 'Z':
  case 'H':
    // The alignment::block adds a sentinel, so this string is NUL-terminated
    // even if the alignment record is malformed.
    return 2 + 1 + strlen(data) + 1;

  default:
    throw bad_format(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' has invalid type ('"
	<< type_ << "')");
  }
}

static string stringize(int) {
  return "FIXME"; // FIXME
}

string alignment::tagfield::value() const {
  switch (type_) {
  case 'A':
    return string(1, data[0]);

  case 'c':
  case 'C':
  case 's':
  case 'S':
  case 'i':
  case 'I':
    return stringize(value_int());

  case 'f':
  case 'd':
    throw std::logic_error("Implement tagfield::value(f/d)"); // TODO

  case 'Z':
  case 'H':
    return string(data);

  default:
    throw bad_format(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' has invalid type ('"
	<< type_ << "')");
  }
}

int alignment::tagfield::value_int() const {
  switch (type_) {
  case 'c':  { signed char   value = data[0]; return value; }
  case 'C':  { unsigned char value = data[0]; return value; }
  case 's':  return convert::int16(data);
  case 'S':  return convert::uint16(data);
  case 'i':  return convert::int32(data);
  // FIXME do something about the big ones...
  case 'I':  return convert::uint32(data);

  case 'f':
  case 'd':
  case 'A':
  case 'Z':
  case 'H':
    // FIXME What exception should we be throwing?
    throw exception(make_string()
	<< "Aux field '" << tag_[0] << tag_[1]
	<< "' is of non-integral type ('" << type_ << "')");

  default:
    throw bad_format(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' has invalid type ('"
	<< type_ << "')");
  }
}

alignment::iterator alignment::find(const char* key) {
  for (iterator it = begin(); it != end(); ++it)
    if (it->tag_equals(key))  return it;

  return end();
}

alignment::const_iterator alignment::find(const char* key) const {
  for (const_iterator it = begin(); it != end(); ++it)
    if (it->tag_equals(key))  return it;

  return end();
}

int alignment::erase(const char* key) {
  int count = 0;
  iterator it = begin();
  while (it != end())
    if (it->tag_equals(key))  it = erase(it), count++;
    else  ++it;

  return count;
}

// Replace the fields in  [start,limit)  by  gap_length  bytes, to be filled in
// as a new auxiliary field by the caller.  Returns  start  (perhaps adjusted);
// if the alignment::block needed to be resized, all other iterators etc are
// invalidated.
// FIXME char* v iterators in doco above
char* alignment::replace_gap(char* start, char* limit, int gap_length) {
  int SENLEN = 0;  // FIXME TODO What about this sentinel?

  int delta = gap_length - (limit - start);
  int newsize = p->size() + delta + SENLEN;

  if (p->capacity() < newsize) {
    char* olddata = p->data();
    resize_unshare_copy(newsize);
    char* newdata = p->data();

    start = &newdata[start - olddata];
    limit = &newdata[limit - olddata];
  }

  memmove(start + gap_length, limit, end().ptr - limit + SENLEN);
  p->c.rest_length += delta;

  return start;
}

alignment::iterator
alignment::replace_string(iterator start, iterator limit,
	    const char* tag, char type, const char* value, int value_length) {
  iterator it = replace_gap(start, limit, 2 + 1 + value_length + 1);

  if (tag)
    it->tag_[0] = tag[0], it->tag_[1] = tag[1];
  it->type_ = type;
  memcpy(it->data, value, value_length);
  it->data[value_length] = '\0';

  return it;
}

alignment::iterator
alignment::replace_(iterator start, iterator limit,
		    const char* tag, int value) {
  iterator it;

  // Pick the shortest representation that can hold the given value,
  // preferring signed to unsigned.

  if (value >= 0) {
    if (value <= INT8_MAX) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(int8_t));
      it->type_ = 'c';
      it->data[0] = value;
    }
    else if (value <= UINT8_MAX) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(uint8_t));
      it->type_ = 'C';
      it->data[0] = value;
    }
    else if (value <= INT16_MAX) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(int16_t));
      it->type_ = 's';
      convert::set_bam_int16(it->data, value);
    }
    else if (value <= UINT16_MAX) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(uint16_t));
      it->type_ = 'S';
      convert::set_bam_uint16(it->data, value);
    }
    else if (value <= INT32_MAX) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(int32_t));
      it->type_ = 'i';
      convert::set_bam_int32(it->data, value);
    }
    else if (value <= int(UINT32_MAX)) { // FIXME signedness of the constant...
      it = replace_gap(start, limit, 2 + 1 + sizeof(uint32_t));
      it->type_ = 'I';
      convert::set_bam_uint32(it->data, value);
    }
    else
      throw std::range_error(make_string() << "Integer aux field '"
	  << tag[0] << tag[1] << "' is out of range");
  }
  else {
    if (value >= INT8_MIN) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(int8_t));
      it->type_ = 'c';
      it->data[0] = value;
    }
    else if (value >= INT16_MIN) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(int16_t));
      it->type_ = 's';
      convert::set_bam_int16(it->data, value);
    }
    else if (value >= INT32_MIN) {
      it = replace_gap(start, limit, 2 + 1 + sizeof(int32_t));
      it->type_ = 'i';
      convert::set_bam_int32(it->data, value);
    }
    else
      throw std::range_error(make_string() << "Integer aux field '"
	  << tag[0] << tag[1] << "' is out of range");
  }

  if (tag)
    it->tag_[0] = tag[0], it->tag_[1] = tag[1];

  return it;
}

alignment::iterator
alignment::replace_(iterator start, iterator limit,
		    const char* tag, const_iterator value) {
  int value_size = value->size();

  iterator it = replace_gap(start, limit, value_size);

  if (tag == NULL)  tag = value->tag_;
  it->tag_[0] = tag[0], it->tag_[1] = tag[1];

  it->type_ = value->type_;
  memcpy(it->data, value->data, value_size - (2 + 1));

  return it;
}

} // namespace sam
