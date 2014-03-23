/*  alignment.cpp -- Classes and functions for SAM/BAM alignment records.

    Copyright (C) 2010-2014 Genome Research Ltd.

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

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#if !(defined UINT8_MAX && defined UINT16_MAX && defined INT32_MAX)
#error - This implementation requires 8-, 16-, and 32-bit integer types
#endif

#include "cansam/sam/alignment.h"

#include <stdexcept>
#include <string>
#include <cstring>
#include <cstdlib>  // for strtol()
#include <climits>  // for UCHAR_MAX

#include <iostream> // FIXME NUKE-ME

#include "cansam/sam/header.h"
#include "cansam/exception.h"
#include "lib/utilities.h"
#include "lib/wire.h"

using std::string;

namespace sam {

// 1. Accessors
//=============

// For order(); 0, FIRST_IN_PAIR, SECOND_IN_PAIR, FIRST_IN_PAIR|SECOND_IN_PAIR.
const int alignment::order_value[4] = { 0, -1, +1, 0 };

string alignment::cigar() const {
  // TODO Probably will build some kind of sam::cigar class...

  int cigar_length = p->c.cigar_length;

  if (cigar_length == 0)
    return "*";

  string text;
  char buffer[format::buffer<int>::size];
  const char* cigar_data = p->cigar_data();
  for (int i = 0; i < cigar_length; i++, cigar_data += sizeof(uint32_t)) {
    uint32_t code = convert::uint32(cigar_data);
    text.append(buffer, format::decimal(buffer, code >> 4) - buffer);
    text += "MIDNSHP=X???????"[code & 0xf];
  }

  return text;
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

namespace {
  const char decode_seq[] =
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
}

// FIXME Use inline template to avoid writing the code twice

char* alignment::unpack_seq(char* dest, const char* raw_seq, int seq_length) {
  const unsigned char* packed = reinterpret_cast<const unsigned char*>(raw_seq);

  int even_length = seq_length & ~1;
  for (int i = 0; i < even_length; i += 2) {
    int ndx = 2 * *packed++;
    *dest++ = decode_seq[ndx];
    *dest++ = decode_seq[ndx+1];
  }

  if (even_length < seq_length) {
    int ndx = 2 * *packed;
    *dest++ = decode_seq[ndx];
  }

  return dest;
}

// Unpack the BAM-encoded sequence data, which is of course not necessarily
// NUL-terminated.  Assumes that DEST is already of sufficient length for
// all this (c.f. unpack_seq()).
void alignment::unpack_seq(string::iterator dest,
			   const char* raw_seq, int seq_length) {
  const unsigned char* packed = reinterpret_cast<const unsigned char*>(raw_seq);

  int even_length = seq_length & ~1;
  for (int i = 0; i < even_length; i += 2) {
    int ndx = 2 * *packed++;
    *dest++ = decode_seq[ndx];
    *dest++ = decode_seq[ndx+1];
  }

  if (even_length < seq_length) {
    int ndx = 2 * *packed;
    *dest = decode_seq[ndx];
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

char* alignment::unpack_qual(char* dest, const char* phred, int seq_length) {
  const char* phredlim = &phred[seq_length];
  while (phred < phredlim) {
    char q = *phred++;
    if (q < 0)  q = 0;
    else if (q > 126 - 33)  q = 126 - 33;
    *dest++ = q + 33;
  }

  return dest;
}

void alignment::unpack_qual(string::iterator dest,
			    const char* phred, int seq_length) {
  for (int i = 0; i < seq_length; i++) {
    char q = *phred++;
    if (q < 0)  q = 0;
    else if (q > 126 - 33)  q = 126 - 33;

    *dest++ = q + 33;
  }
}

// FIXME Hmmm... static method of cigar class?
int cigar_operator_count(const char* s) {
  int n = 0;

  if (s[0] == '*' && s[1] == '\0') {
    // An unspecified cigar string, written as "*", is represented as
    // empty in BAM, so leave  n  as 0.
  }
  else {
    while (*s)
      if (! isdigit(*s++))  n++;
  }

  return n;
}

// FIXME Either should be static alignment::pack_cigar(), or hello cigar class
// FIXME Think about how much this can write and whether it throws on all cigar
// syntax errors; and whether cigar_operator_count() is named right and is
// always a bound on how much this can write.
void pack_cigar(char* dest, const char* cigar) {
  static const string cigar_operators = "MIDNSHP=X";

  const char *s = cigar;
  while (*s) {
    const char* sdigit = s;
    int len = 0;
    while (isdigit(*s))
      len = 10 * len + *s++ - '0';

    if (s == sdigit)
      throw bad_format(make_string()
	  << "Missing digits in CIGAR string ('" << cigar << "')");

    char op_char = *s++;
    size_t op = cigar_operators.find(op_char);
    if (op == string::npos) {
      if (op_char == '\0')
	throw bad_format(make_string()
	    << "Truncated CIGAR string ('" << cigar << "')");
      else
	throw bad_format(make_string()
	    << "Invalid operator ('" << op_char << "') in CIGAR string ('"
	    << cigar << "')");
    }

    convert::set_bam_uint32(dest, (len << 4) | op);
    dest += sizeof(uint32_t);
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
  char* qbuffer = replace_gap(p->name_data(),
		      p->name_data() + p->c.name_length, qname_length + 1);

  p->c.name_length = qname_length + 1;
  memcpy(qbuffer, qname_data, qname_length);
  qbuffer[qname_length] = '\0';
}

void alignment::set_flags(int flags) {
  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.flags = flags;
}

void alignment::set_rindex(int rindex) {
  int ref_size = collection::find(p->h.cindex).ref_size();
  if (rindex < -1 || rindex >= ref_size)
    throw std::out_of_range(make_string()
	<< "New rindex value (" << rindex
	<< ") is outwith range [-1," << ref_size << ")");

  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.rindex = rindex;
}

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

void alignment::set_cigar(const char* cigar) {
  int new_cigar_length = cigar_operator_count(cigar);
  char* cbuffer = replace_gap(p->cigar_data(),
		      p->cigar_data() + sizeof(uint32_t) * p->c.cigar_length,
		      sizeof(uint32_t) * new_cigar_length);

  p->c.cigar_length = new_cigar_length;
  if (! (cigar[0] == '*' && cigar[1] == '\0'))
    pack_cigar(cbuffer, cigar);
}

void alignment::set_mate_rindex(int rindex) {
  int ref_size = collection::find(p->h.cindex).ref_size();
  if (rindex < -1 || rindex >= ref_size)
    throw std::out_of_range(make_string()
	<< "New mate_rindex value (" << rindex
	<< ") is outwith range [-1," << ref_size << ")");

  if (p == &empty_block)  resize_unshare_copy(p->size());
  p->c.mate_rindex = rindex;
}

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
alignment::block alignment::empty_block = {
  { 0 /* 37, if truth be told */, 0 },
  { 33, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0 },
  { '\0' /* an empty qname C-string */ }
};

// Allocate a new block and initialise its capacity field.
alignment::block* alignment::block::create(int payload_size) {
  char* cp = new char[sizeof(block_header) + payload_size];
  block* p = reinterpret_cast<block*>(cp);

  p->h.capacity = payload_size;
  return p;
}

// Copy the block contents (but don't overwrite DEST's capacity field),
// assuming that the destination's capacity suffices for the source's size.
void alignment::block::copy(block* dest, const block* src) {
  memcpy(&dest->h.cindex, &src->h.cindex, sizeof(src->h.cindex) + src->size());
}

// Deallocate the block (which must not be empty_block).
void alignment::block::destroy(block* p) {
  char* cp = reinterpret_cast<char*>(p);
  delete [] cp;
}

/* Resize the alignment's block, unsharing it if it is currently the shared
empty block, while maintaining the existing contents.  Used by mutators to
ensure the block is sufficiently sized and writable.  */
void alignment::resize_unshare_copy(int payload_size) {
  block* newp = block::create(payload_size);
  block* oldp = p;
  block::copy(newp, oldp);
  p = newp;
  if (oldp != &empty_block)  block::destroy(oldp);
}

/* Resize the alignment's block, unsharing it if it is currently the shared
empty block, but not maintaining the existing contents.  Used by assignments
to ensure the block is sufficiently sized and writable.  */
void alignment::resize_unshare_discard(int payload_size) {
  block* newp = block::create(payload_size);
  block* oldp = p;
  p = newp;
  if (oldp != &empty_block)  block::destroy(oldp);
}

/* This assignment operator simply copies the pointed-to block, but must first
ensure that the destination is large enough and is not the shared empty block,
and must explicitly check for the self-assignment case.  */
// FIXME Should we de-pessimize the empty = empty case?
alignment& alignment::operator= (const alignment& aln) {
  int newsize = aln.p->size();
  if (p->capacity() < newsize)
    resize_unshare_discard(newsize);
  else if (p == aln.p)
    return *this;

  block::copy(p, aln.p);
  return *this;
}

alignment::alignment(const alignment& aln) {
  if (aln.p == &empty_block)
    p = &empty_block;
  else {
    p = block::create(aln.p->size());
    block::copy(p, aln.p);
  }
}

#if 1
alignment& alignment::assign(const string& /*line*/) {
  // FIXME
  return *this;
}
#endif

// FIXME This is pretty crap, and ought to be in utilities.h
int decimal(const char* s, const char* field_name) {
  int x;
  if (*parse::decimal(s, x) != '\0')
    throw bad_format(make_string()
	<< "Trailing gunge in decimal " << field_name
	<< " field ('" << s << "')");
  return x;
}

int parse_flags(const char* s) {
  int value;

  if (isdigit(s[0])) {
    char* slim;
    value = strtol(s, &slim, 0);
    if (*slim != '\0')
      throw bad_format(make_string()
	  << "Flag value ('" << s << "') is non-numeric");
  }
  else {
    value = 0;
    for (const char* schr = s; *schr; schr++)
      switch (*schr) {
      case 'p':  value |= PAIRED;  break;
      case 'P':  value |= PROPER_PAIRED;  break;
      case 'u':  value |= UNMAPPED;  break;
      case 'U':  value |= MATE_UNMAPPED;  break;
      case 'r':  value |= REVERSE_STRAND;  break;
      case 'R':  value |= MATE_REVERSE_STRAND;  break;
      case '1':  value |= FIRST_IN_PAIR;  break;
      case '2':  value |= SECOND_IN_PAIR;  break;
      case 's':  value |= NONPRIMARY;  break;
      case 'q':  value |= QUALITY_FAILED;  break;
      case 'd':  value |= DUPLICATE;  break;
      case 'S':  value |= SUPPLEMENTARY;  break;

      case 'f':
      case 'F':
      case '_':
	break;

      default:
	throw bad_format(make_string()
	    << "Unknown flag character ('" << *schr << "') in value ('" << s
	    << "')");
      }
  }

  return value;
}

void parse_flags(const std::string& text, int& positive, int& negative) {
  size_t next_pos;
  for (size_t pos = 0; pos < text.length(); pos = next_pos) {
    next_pos = text.find_first_of("+-", pos + 1);

    bool positive_item = (text[pos] != '-');
    if (text[pos] == '+' || text[pos] == '-')  pos++;
    string item = text.substr(pos, next_pos - pos);

    int& flags  = positive_item? positive : negative;
    int& others = positive_item? negative : positive;

    flags |= parse_flags(item);
    if (item.find('f') != string::npos)  others |= REVERSE_STRAND;
    if (item.find('F') != string::npos)  others |= MATE_REVERSE_STRAND;
  }
}

int alignment::sam_length() const {
  int len = 0;

  len += p->c.name_length - 1;  // name_length includes the trailing NUL

  if (flags() & (PAIRED | PROPER_PAIRED | MATE_UNMAPPED | MATE_REVERSE_STRAND |
		 FIRST_IN_PAIR | SECOND_IN_PAIR)) {
    len += 1 + 12;  // Longest FLAG is "urURpP12sSqd"
  }
  else
    len += 1 + 6;  // Longest unpaired FLAG is "ursSqd"

  len += 1 + rname().length();
  len += 1 + format::buffer<coord_t>::size;
  len += 1 + 3;  // Longest MAPQ is "255"
  // TODO Look at a reasonable bound on the numbers in the CIGAR string
  len += 1 + p->c.cigar_length * (format::buffer<int>::size + 1);
  len += 1 + mate_rname().length();
  len += 1 + format::buffer<coord_t>::size;
  len += 1 + format::buffer<scoord_t>::size;
  len += 1 + p->c.read_length;  // SEQ
  len += 1 + p->c.read_length;  // QUAL

  for (const_iterator it = begin(); it != end(); ++it)
    len += 1 + it->sam_length();

  return len;
}

// As strcpy(), but returns the first unused position in DEST.
inline char* copy(char* dest, const char* s) {
  char c;
  while ((c = *s++) != '\0')  *dest++ = c;
  return dest;
}

char* format_sam(char* dest, const alignment& aln, const std::ios& format) {
  dest = copy(dest, aln.qname_c_str());

  *dest++ = '\t';
  dest = format_sam(dest, aln.flags(), format);

  *dest++ = '\t';
  if (aln.rindex() < 0)  *dest++ = '*';
  else  dest = copy(dest, aln.rname_c_str());

  *dest++ = '\t';
  dest = format::decimal(dest, aln.pos());

  *dest++ = '\t';
  dest = format::decimal(dest, aln.mapq());

  *dest++ = '\t';
  if (aln.p->c.cigar_length == 0)
    *dest++ = '*';
  else {
    // TODO This too probably will end up in some kind of sam::cigar class
    const char* cigar;
    const char* cigarlim = aln.p->cigar_data() +
			      sizeof(uint32_t) * aln.p->c.cigar_length;
    for (cigar = aln.p->cigar_data(); cigar < cigarlim;
	 cigar += sizeof(uint32_t)) {
      uint32_t code = convert::uint32(cigar);
      dest = format::decimal(dest, code >> 4);
      *dest++ = "MIDNSHP=X???????"[code & 0xf];
    }
  }

  *dest++ = '\t';
  if (aln.mate_rindex() < 0)  *dest++ = '*';
  else if (aln.mate_rindex() == aln.rindex())  *dest++ = '=';
  else  dest = copy(dest, aln.mate_rname_c_str());

  *dest++ = '\t';
  // FIXME Do we need to do anything special for 0 or -1 i.e. unmapped?
  dest = format::decimal(dest, aln.mate_pos());

  *dest++ = '\t';
  dest = format::decimal(dest, aln.isize());

  *dest++ = '\t';
  if (aln.length() == 0)  *dest++ = '*';
  else  dest = alignment::unpack_seq(dest, aln.seq_raw_data(), aln.length());

  *dest++ = '\t';
  if (aln.length() == 0 || aln.qual_raw_data()[0] == '\xff')  *dest++ = '*';
  else  dest = alignment::unpack_qual(dest, aln.qual_raw_data(), aln.length());

  for (alignment::const_iterator it = aln.begin(); it != aln.end(); ++it) {
    *dest++ = '\t';
    dest = format_sam(dest, *it);
  }

  return dest;
}

char* format_sam(char* dest, int flags, const std::ios& format) {
  std::ios::fmtflags fmtflags = format.flags();

  if (fmtflags & std::ios::boolalpha) {
    if (flags & UNMAPPED)  *dest++ = 'u';
    if (flags & REVERSE_STRAND)  *dest++ = 'r';
    else if (! (flags & UNMAPPED))  *dest++ = 'f';

    if (flags & MATE_UNMAPPED)  *dest++ = 'U';
    if (flags & MATE_REVERSE_STRAND)  *dest++ = 'R';
    else if ((flags & (PAIRED | MATE_UNMAPPED)) == PAIRED)  *dest++ = 'F';

    if (flags & PAIRED)          *dest++ = 'p';
    if (flags & PROPER_PAIRED)   *dest++ = 'P';
    if (flags & FIRST_IN_PAIR)   *dest++ = '1';
    if (flags & SECOND_IN_PAIR)  *dest++ = '2';
    if (flags & NONPRIMARY)      *dest++ = 's';
    if (flags & SUPPLEMENTARY)   *dest++ = 'S';
    if (flags & QUALITY_FAILED)  *dest++ = 'q';
    if (flags & DUPLICATE)       *dest++ = 'd';
  }
  else if (fmtflags & std::ios::oct)
    dest = format::octal(dest, flags);
  else if (fmtflags & std::ios::hex)
    dest = format::hexadecimal(dest, flags);
  else
    dest = format::decimal(dest, flags);

  return dest;
}

char* format_sam(char* dest, const alignment::tagfield& aux) {
  *dest++ = aux.tag_[0], *dest++ = aux.tag_[1];
  *dest++ = ':';
  switch (aux.type_) {
  case 'A':
    *dest++ = 'A';
    *dest++ = ':';
    *dest++ = aux.data[0];
    break;

  case 'c':
  case 's':
  case 'i':
    *dest++ = 'i';
    *dest++ = ':';
    dest = format::decimal(dest, aux.value<int>());
    break;

  case 'C':
  case 'S':
  case 'I':
    // FIXME These ones should be value_uint() or so (especially I!)
    *dest++ = 'i';
    *dest++ = ':';
    dest = format::decimal(dest, aux.value<int>());
    break;

  case 'f':
  case 'd':
    throw std::logic_error("Implement formatsam() for tagfield(f/d)"); // TODO

  case 'Z':
  case 'H':
    *dest++ = aux.type_;
    *dest++ = ':';
    dest = copy(dest, aux.data);
    break;

  default:
    throw bad_format(make_string()
	<< "Aux field '" << aux.tag_[0] << aux.tag_[1]
	<< "' has invalid type ('" << aux.type_ << "')");
  }

  return dest;
}

void alignment::assign(int nfields, const std::vector<char*>& fields, int cindex) {
  // An alignment in SAM format is a tab-separated line containing fields
  // ordered as:
  // qname flag rname pos mapq cigar mrname mpos isize seq qual aux...
  // 0     1    2     3   4    5     6      7    8     9   10   11...
  enum { qname, flag, rname, pos, mapq, cigar, mrname, mpos, isize, seq, qual,
	 firstaux };

  if (nfields < 11)
    throw bad_format("Too few fields in SAM record");

  int name_length = fields[qname+1] - fields[qname]; // including terminator
  int cigar_length = cigar_operator_count(fields[cigar]);
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

  int size = sizeof(bamcore) + name_length + (cigar_length * sizeof(uint32_t))
	     + (seq_length+1)/2 + seq_length;

  int aux_size = 0;
  // These aux fields will be properly validated in push_back_sam() below.
  for (int i = firstaux; i < nfields; i++)
    aux_size +=
	alignment::tagfield::size_sam(fields[i], fields[i+1] - fields[i] - 1);

  if (p->capacity() < size + aux_size)
    resize_unshare_discard(size + aux_size);

  p->h.cindex = cindex;
  collection& collection = collection::find(cindex);

  // The aux space will be added to rest_length during push_back_sam() below.
  p->c.rest_length = size - sizeof(p->c.rest_length);

  p->c.rindex = collection.findseq(fields[rname]).index(); // a name or "*"
  p->c.zpos = decimal(fields[pos], "POS") - 1; // a 1-based pos or 0
  p->c.name_length = name_length;
  p->c.mapq = decimal(fields[mapq], "MAPQ"); // 0..255
  p->c.bin = -1;
  p->c.cigar_length = cigar_length;
  p->c.flags = parse_flags(fields[flag]);
  p->c.read_length = seq_length;

  //a name or "*" or "="
  if (fields[mrname][0] == '=' && fields[mrname][1] == '\0')
    p->c.mate_rindex = p->c.rindex;
  else
    p->c.mate_rindex = collection.findseq(fields[mrname]).index();

  p->c.mate_zpos = decimal(fields[mpos], "MPOS") - 1; // a 1-based pos or 0
  p->c.isize = decimal(fields[isize], "ISIZE"); // an insert size or 0

  memcpy(p->name_data(), fields[qname], name_length);

  if (! (fields[cigar][0] == '*' || fields[cigar][1] == '\0'))
    pack_cigar(p->cigar_data(), fields[cigar]);

  // In BAM: int32_t read_len is given, seq is (read_len+1)/2 bytes;
  //         qual is read_len bytes (Phred qualities), maybe 0xFF x read_len.

  // In SAM: either seq is "*" so read_len is 0, or seq is read_len base chars;
  //         qual is "*" (so 0xFF x read_len) or it's P-33 x read_len qualchars.

  pack_seq(p->seq_data(), fields[seq], seq_length);
  if (qual_length > 0)
    pack_qual(p->qual_data(), fields[qual], seq_length);
  else
    memset(p->qual_data(), '\xff', seq_length);

  for (int i = firstaux; i < nfields; i++)
    push_back_sam(fields[i], fields[i+1] - fields[i] - 1);
}

void alignment::push_back_sam(const char* aux, int length) {
  // A SAM-formatted aux field is in the format "TG:T:[VALUE]".

  if (length < 5 || aux[2] != ':' || aux[4] != ':')
    throw bad_format("Malformatted aux field");

//std::clog << "push_back_sam(\"" << aux << "\", " << length << ")\n";
  const char* tag = aux;
  const char* value = &aux[5];
  int value_length = length - 5;

  switch (aux[3]) {
  case 'A':
    if (value_length != 1)
      throw bad_format("Type 'A' aux field has length other than 1");
    push_back(tag, value[0]);
    break;

  case 'i': {
    // TODO  Conceivably take heroic steps to represent unsigned [2^31, 2^32)
    int ivalue;
    if (*parse::decimal(value, ivalue) != '\0')
      throw bad_format(make_string()
	  << "Numeric aux field has non-numeric value '" << value << "')");
    push_back(tag, ivalue);
    }
    break;

  case 'f':
    throw std::logic_error("Aux 'f' field not implemented");  // TODO

  case 'd':
    throw std::logic_error("Aux 'd' field not implemented");  // TODO

  case 'Z':
    replace_string(end(), end(), tag, 'Z', value, value_length);
    break;

  case 'H':
    throw std::logic_error("Aux 'H' field not implemented");  // TODO

  default:
    throw bad_format(make_string()
	<< "Aux field '" << tag[0] << tag[1] << "' has invalid type ('"
	<< aux[3] << "') for SAM format");
  }
}


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

int alignment::tagfield::sam_length() const {
  // Returns the number of characters in "TG:T:VALUE", i.e., 5 + VALUE_length.
  // For integer types, this is a conservative (i.e., generous) approximation.
  // As this is primarily used for sizing output buffers, there is no need to
  // waste time looking at the value to determine a more exact length.

  switch (type_) {
  case 'A':  return 5 + 1;
  case 'c':  return 5 + 4;   // "-128" and similar are the longest possible
  case 'C':  return 5 + 3;   // "255"
  case 's':  return 5 + 6;   // "-32768"
  case 'S':  return 5 + 5;   // "65535"
  case 'i':  return 5 + 11;  // "-2147483648"
  case 'I':  return 5 + 10;  // "4294967296"

  case 'f':  throw std::logic_error("Aux 'f' field not implemented");  // TODO
  case 'd':  throw std::logic_error("Aux 'd' field not implemented");  // TODO

  case 'Z':
  case 'H':
    return 5 + strlen(data);

  default:
    throw bad_format(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' has invalid type ('"
	<< type_ << "')");
  }
}

int alignment::tagfield::size_sam(const char* text, int text_length) {
  // Return 0 (invalid) if the text is clearly not of the form "TG:T:[VALUE]".
  if (text_length < 5)  return 0;

  const char* value = &text[5];
  int length = text_length - 5;

  switch (text[3]) {
  case 'A':
    return 2 + 1 + 1;

  case 'i':
    if (value[0] != '-') {
      if (length <= 2)  // <= 99
	return 2 + 1 + sizeof(uint8_t);
      else if (length <= 4 || (length == 5 && value[0] <= '5'))  // <= 59999
	return 2 + 1 + sizeof(uint16_t);
      else
	return 2 + 1 + sizeof(uint32_t);
    }
    else {
      if (length <= 3)  // >= -99
	return 2 + 1 + sizeof(int8_t);
      else if (length <= 5 || (length == 6 && value[1] <= '2'))  // >= -29999
	return 2 + 1 + sizeof(int16_t);
      else
	return 2 + 1 + sizeof(int32_t);
    }

  case 'f':
    return 2 + 1 + 4;

  case 'd':
    return 2 + 1 + 8;

  case 'Z':
  case 'H':
    return 2 + 1 + length + 1;

  default:
    return 0;
  }
}

string& alignment::tagfield::value(string& dest) const {
  switch (type_) {
  case 'A':
    dest = data[0];
    break;

  case 'c':
  case 's':
  case 'i': {
    char buffer[format::buffer<int>::size];
    dest.assign(buffer, format::decimal(buffer, value<int>()) - buffer);
    break;
    }

  case 'C':
  case 'S':
  case 'I': {
    // FIXME These ones should be value_uint() or so (especially I!)
    char buffer[format::buffer<int>::size];
    dest.assign(buffer, format::decimal(buffer, value<int>()) - buffer);
    break;
    }

  case 'f':
  case 'd':
    throw std::logic_error("Implement tagfield::value(f/d)"); // TODO

  case 'Z':
  case 'H':
    dest = data;
    break;

  default:
    throw bad_format(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' has invalid type ('"
	<< type_ << "')");
  }

  return dest;
}

template<> const char* alignment::tagfield::value() const {
  switch (type_) {
  case 'Z':
  case 'H':
    return data;

  case 'c':  case 'C':
  case 's':  case 'S':
  case 'i':  case 'I':
  case 'f':  case 'd':
  case 'A':
    throw sam::exception(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' is of non-string type ('"
	<< type_ << "')");

  default:
    throw bad_format(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' has invalid type ('"
	<< type_ << "')");
  }
}

template<> int alignment::tagfield::value() const {
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

template<> char alignment::tagfield::value() const {
  switch (type_) {
  case 'A':
    return data[0];

  case 'Z':
    if (data[0] != '\0' && data[1] == '\0')
      return data[0];
    else
      throw sam::exception(make_string()
	  << "Aux field '" << tag_[0] << tag_[1]
	  << "' has length other than 1");

  case 'c':  case 'C':
  case 's':  case 'S':
  case 'i':  case 'I':
  case 'f':  case 'd':
  case 'H':
    throw sam::exception(make_string()
	<< "Aux field '" << tag_[0] << tag_[1] << "' is of non-char type ('"
	<< type_ << "')");

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

alignment::const_iterator alignment::find_or_throw(const char* key) const {
  for (const_iterator it = begin(); it != end(); ++it)
    if (it->tag_equals(key))  return it;

  throw sam::exception(make_string()
      << "Aux field '" << key[0] << key[1] << "' not found");
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
		    const char* tag, char value) {
  iterator it = replace_gap(start, limit, 2 + 1 + 1);

  if (tag)
    it->tag_[0] = tag[0], it->tag_[1] = tag[1];
  it->type_ = 'A';
  it->data[0] = value;

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
