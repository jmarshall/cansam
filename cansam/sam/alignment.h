/// @file cansam/sam/alignment.h
/// Classes and functions for SAM/BAM alignment records

/*  Copyright (C) 2010-2014 Genome Research Ltd.
    Portions copyright (C) 2020 University of Glasgow.

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

#ifndef CANSAM_SAM_ALIGNMENT_H
#define CANSAM_SAM_ALIGNMENT_H

#include <string>
#include <vector>
#include <iterator>
//#include <algorithm>  // FIXME NUKE-ME, but for specialising std::swap IIRC
#include <iosfwd>
#include <cstring>

#include <stdint.h>

#include "cansam/types.h"
#include "cansam/sam/header.h"

namespace sam {

class bamio;
class samio;

// FIXME do something about namespaces so you don't have to say sam::PAIRED etc

/** The @e flags field of an alignment contains a number of flag bits pertaining
    to the read that has been aligned and a few pertaining to its mate.  */
enum alignment_flag {
  PAIRED              = 0x001, ///< The read is paired with a mate
  PROPER_PAIRED       = 0x002, ///< The read and its mate are mapped as a proper pair
  UNMAPPED            = 0x004, ///< The read is unmapped
  MATE_UNMAPPED       = 0x008, ///< The mate is unmapped
  REVERSE_STRAND      = 0x010, ///< The read's strand (0: forward; 1: reverse)
  MATE_REVERSE_STRAND = 0x020, ///< The mate's strand (0: forward; 1: reverse)
  FIRST_IN_PAIR       = 0x040, ///< The read is the first of a pair
  SECOND_IN_PAIR      = 0x080, ///< The read is the second of a pair
  NONPRIMARY          = 0x100, ///< This alignment is not the primary one for the read
  QUALITY_FAILED      = 0x200, ///< Failed platform or vendor quality checks
  DUPLICATE           = 0x400, ///< PCR duplicate or optical duplicate
  SUPPLEMENTARY       = 0x800  ///< This alignment is a supplementary one for the read
};

/// Returns the bitwise representation of @a flags
/** @param flags A string representing a set of flags either numerically
or symbolically.
@return The bitwise integer representation of the alignment flags specified
by @a flags, which may be either decimal, octal, hexadecimal, or symbolic.
@relatesalso alignment */
int parse_flags(const char* flags);

/// Returns the bitwise representation of @a flags
/** @relatesalso alignment */
inline int parse_flags(const std::string& flags)
  { return parse_flags(flags.c_str()); }

/// Accumulates @a flags into signed categories
/** @param flags A string representing several sets of alignment flags,
separated by '+' or '-' characters.
@param positive Gains set bits corresponding to flags following '+' characters.
@param negative Gains set bits corresponding to flags following '-' characters.
@relatesalso alignment */
void parse_flags(const std::string& flags, int& positive, int& negative);

/// Returns the BAM bin number (1-based)
/** Returns the BAM bin number for an alignment spanning [@a pos, @a right],
    i.e., a 1-based range.  */
int calc_bin(coord_t pos, coord_t right);

/// Returns the BAM bin number (0-based)
/** Returns the BAM bin number for an alignment spanning [@a zpos, @a zright],
    i.e., a 0-based range.  */
int calc_zbin(coord_t zpos, coord_t zright);

/** CIGAR string operators.  */
enum cigar_opcode {
  MATCH,          ///< Alignment match (can be a sequence match or mismatch)
  INSERTION,      ///< Insertion to the reference
  DELETION,       ///< Deletion from the reference
  REF_SKIP,       ///< Skipped region from the reference
  SOFT_CLIP,      ///< Soft clipping (clipped sequences present in SEQ)
  HARD_CLIP,      ///< Hard clipping (clipped sequences NOT present in SEQ)
  PADDING,        ///< Padding (silent deletion from padded reference)
  MATCH_EQUAL,    ///< Sequence match
  MATCH_DIFF      ///< Sequence mismatch
};

/** @class sam::cigar_op cansam/sam/alignment.h
*/
class cigar_op {
public:
  cigar_op(int length, char opchar) : data_(length << 4 | encode(opchar)) { }
  cigar_op(const cigar_op& cigar) : data_(cigar.data_) { }
  ~cigar_op() { }

  cigar_op& operator= (const cigar_op& cigar)
    { data_ = cigar.data_; return *this; }

  cigar_opcode opcode() const { return static_cast<cigar_opcode>(data_ & 0xf); }
  char opchar() const { return opchars[data_ & 0xf]; }
  int length() const { return data_ >> 4; }

  static cigar_opcode encode(char opchar);
  static char decode(cigar_opcode opcode) { return opchars[opcode]; }

private:
  // @cond private
  friend class alignment;
  static const char opchars[16];

  cigar_op(const char* ptr);

  uint32_t data_;
  // @endcond
};

/// Print a CIGAR operation to the stream
/** Writes a CIGAR record to a stream as text in SAM format.
@relatesalso cigar_op */
std::ostream& operator<< (std::ostream& stream, const cigar_op& cigar);

/// Write a CIGAR operation to @a dest in SAM format
/** @param dest  Character array to be written to; must have space for
at least 10 characters.
@param cigar  The CIGAR operation to be formatted.
@return  A pointer to the first unused character position in @a dest.
@relatesalso cigar_op */
char* format_sam(char* dest, const cigar_op& cigar);

/// Print an unpacked CIGAR string to the stream
/** Writes a vector of CIGAR records to a stream as text in SAM format.
@relatesalso cigar_op */
std::ostream& operator<< (std::ostream& stream, const std::vector<cigar_op>&);

/** @class sam::alignment cansam/sam/alignment.h
    @brief SAM/BAM alignment record

Blah blah blah about alignment records, including how they are represented and
the performance implications for user code, e.g., use swap() in preference to
assignment, and about thread considerations.

The various auxiliary field function templates take a @em ValueType parameter,
which may be of any of the following types:
  - <tt>const std::string&</tt> or <tt>const char*</tt>
  - @c int
  - @c char (eventually)
  - @c float (eventually)
  - @c double (eventually)
  - @c std::vector<uint8_t> (eventually)
  - @c const_iterator (within a different alignment object)

Note that a @c const_iterator used as a source value must be a dereferenceable
iterator pointing to an alignment::tagfield in a @em different alignment object
to the one being modified.
*/
class alignment {
public:
  /// Construct an empty alignment
  alignment() : p(&empty_block) { }

#if 0
  /// Construct an alignment by splitting up a tab-separated text string
  // FIXME what reference_thingie does this use?
  explicit alignment(const std::string& line) { assign(line); }
#endif

  /// Construct a copy of an alignment
  alignment(const alignment& aln);

  //  Destroy this alignment object (not interesting enough to warrant ///)
  ~alignment() { if (p != &empty_block)  block::destroy(p); }

  /// Copy an alignment
  alignment& operator= (const alignment& aln);

  /// Assign to this alignment by splitting up a tab-separated text string
  // FIXME what reference_thingie does this use?
  alignment& assign(const std::string& line);

  /// Swap this alignment with another
  void swap(alignment& aln) { block* tmp = p; p = aln.p; aln.p = tmp; }

#if 0
  /// Ensure this alignment has space for a BAM-formatted record
  /// of @a size bytes
  // FIXME Surely this can be worded better...
  void reserve(int size)
    { if (p->capacity() < size)  resize_unshare_copy(size); }
#endif

  /// Returns the approximate number of characters in the SAM representation
  /// of this alignment
  int sam_length() const;

  /** @name Field accessors
  Two variants are provided for the @em POS and @em MPOS fields: @c %pos()
  et al return 1-based coordinates, while @c %zpos() et al return the same
  positions in 0-based coordinates.  */
  //@{
  std::string qname() const
    { return std::string(p->name_data(), p->c.name_length - 1); }

  int flags() const { return p->c.flags; }

  /** Returns an index uniquely identifying the read's reference sequence
  within the collection of which it is a part.
  @return The index, or -1 if @em RNAME indicates that the read is unmapped. */
  int rindex() const { return p->c.rindex; }

  /** Returns the name of the read's reference sequence.
  @return The name, or "*" if @em RNAME indicates that the read is unmapped. */
  std::string rname() const
    { return collection::find(p->h.cindex).findseq(p->c.rindex).name(); }

  coord_t pos() const  { return p->c.zpos+1; }
  coord_t zpos() const { return p->c.zpos; }
  int mapq() const { return p->c.mapq; }

  /// Number of CIGAR operations
  size_t cigar_length() const { return p->c.cigar_length; }

  cigar_op cigar(size_t i) const;
  template <typename CigarType> CigarType cigar() const;
  std::string& cigar(std::string& dest) const;
  std::vector<cigar_op>& cigar(std::vector<cigar_op>& dest) const;

  /** Returns, for paired reads, an index uniquely identifying the mate read's
  reference sequence within the collection of which it is a part.
  @return The index, or -1 if the read is unpaired or @em MRNM indicates that
  the mate is unmapped.  */
  int mate_rindex() const { return p->c.mate_rindex; }

  /** Returns, for paired reads, the name of the mate read's reference sequence.
  @return The name, or "*" if the read is unpaired or @em MRNM indicates that
  the mate is unmapped.
  @note The actual name is returned, even when this field would appear as "="
  in a SAM file.  */
  std::string mate_rname() const
    { return collection::find(p->h.cindex).findseq(p->c.mate_rindex).name(); }

  coord_t mate_pos() const  { return p->c.mate_zpos+1; }
  coord_t mate_zpos() const { return p->c.mate_zpos; }
  scoord_t isize() const { return p->c.isize; }

  std::string seq() const
    { std::string dest(length(), 'X');
      unpack_seq(dest.begin(), seq_raw_data(), length());
      return dest; }

  std::string qual() const
    { std::string dest(length(), 'X');
      unpack_qual(dest.begin(), qual_raw_data(), length());
      return dest; }

  /// Read length
  int length() const { return p->c.read_length; }

  /// BAM bin number (derived, if necessary, from @em POS and @em CIGAR)
  // FIXME How does this work for unmapped?
  int bin() const
    { if (p->c.bin == unknown_bin)  p->c.bin = calc_zbin(zpos(), right_zpos());
      return p->c.bin; }

  /** Returns the value of the auxiliary field with the given @a tag,
  or throws a sam::exception if the field cannot be expressed as
  a @a ValueType or if there is no such field.
  @note In most cases, it will be necessary to use explicit syntax such
  as <tt>aux<int>("NM")</tt> to specify which type is desired.  */
  template <typename ValueType>
  ValueType aux(const char* tag) const
    { return find_or_throw(tag)->template value<ValueType>(); }

  /** Returns the value of the auxiliary field with the given @a tag,
  or @a default_value if there is no such field; or throws a sam::exception
  if the field is present but cannot be expressed as a @a ValueType.  */
  template <typename ValueType>
  ValueType aux(const char* tag, ValueType default_value) const
    { const_iterator it = find(tag);
      return (it != end())? it->template value<ValueType>() : default_value; }
  //@}

  /** @name Additional field accessors
  Some accessors returning strings also have corresponding versions returning
  C strings or (non-NUL-terminated) @c char arrays.  These versions are cheaper
  as a suitable representation already exists within the alignment object, so
  no @c std::string construction or memory allocation is necessary.

  Pointers returned by an alignment object become invalid when any of
  that object's non-const member functions are subsequently called.  */
  //@{
  /// Query name
  const char* qname_c_str() const { return p->name_data(); }

  /// Query name length (not including the NUL terminator)
  int qname_length() const { return p->c.name_length - 1; }

  /// Assigns query name to @a dest (and returns @a dest)
  std::string& qname(std::string& dest) const
    { return dest.assign(p->name_data(), p->c.name_length - 1); }

  /// Reference name (or @c NULL)
  const char* rname_c_str() const
    { return collection::find(p->h.cindex).findseq(p->c.rindex).name_c_str(); }

  /// Mate's reference name (or @c NULL)
  const char* mate_rname_c_str() const
    { return
	collection::find(p->h.cindex).findseq(p->c.mate_rindex).name_c_str(); }

  /// Assigns sequence to @a dest (and returns @a dest)
  std::string& seq(std::string& dest) const
    { unpack_seq(dest, seq_raw_data(), length()); return dest; }

  /// Sequence (packed as two bases per byte; not NUL-terminated)
  const char* seq_raw_data() const { return p->seq_data(); }

  /// Assigns ASCII quality string to @a dest (and returns @a dest)
  std::string& qual(std::string& dest) const
    { unpack_qual(dest, qual_raw_data(), length()); return dest; }

  /// Quality string (BLAH raw phred scores; not NUL-terminated)
  const char* qual_raw_data() const { return p->qual_data(); }

  std::string& aux(std::string& dest, const char* tag) const
    { return find_or_throw(tag)->value(dest); }

  std::string& aux(std::string& dest,
		   const char* tag, const char* default_value) const
    { const_iterator it = find(tag);
      return (it != end())? it->value(dest) : dest.assign(default_value); }

  std::string& aux(std::string& dest,
		   const char* tag, const std::string& default_value) const
    { const_iterator it = find(tag);
      return (it != end())? it->value(dest) : dest.assign(default_value); }
  //@}

  /** @name Container functionality
  Alignment records provide limited container-style access to their
  auxiliary fields.

  The @c sam::alignment::iterator and @c sam::alignment::const_iterator classes
  are <b>forward iterators</b> providing all the usual iterator functionality:
  copying, assignment, pre- and post-increment, equality and inequality tests,
  and dereferencing, which produces an alignment::tagfield through which the
  pointed-to auxiliary field's properties can be accessed (but not modified;
  see also the iterator variant of set_aux() below).

  The possible types for the @c ValueType parameter are listed below.  */
  //@{
  /** @class sam::alignment::tagfield cansam/sam/alignment.h
      @brief Helper class representing an auxiliary field as seen via an
	     iterator

  Dereferencing a sam::alignment @c iterator or @c const_iterator produces
  (a reference to) an instance of this class.

  @note There are no mutator methods, even if it is a mutable @c iterator that
  has been dereferenced; use sam::alignment::set_aux() to change the value of
  an auxiliary field via an @c iterator.  */
  class tagfield {
  public:
    /// Two-character field tag
    std::string tag() const { return std::string(tag_, sizeof tag_); }

    /// Field type, as it would appear in a BAM file
    /** @return One of 'A', 'Z', 'i', 'c', 'S', etc, including the subtypes
    that only appear in a BAM file.  */
    char type() const { return type_; }

    /// Field value
    template <typename ValueType>
    ValueType value() const;

#if 0
    /// Field value, as it would appear in a SAM file
    /** @return The text string representation of the field's value (rather
    than any binary representation that might appear in a BAM file).  */
    std::string value_str() const;
#endif

    /// Assigns SAM-style field value to @a dest (and returns @a dest)
    std::string& value(std::string& dest) const;

    /// Returns whether the field's tag is the same as the given @a key_tag
    bool tag_equals(const char* key_tag) const
      { return tag_[0] == key_tag[0] && tag_[1] == key_tag[1]; }

    /// Number of bytes in the BAM representation of this field
    int size() const;

    /// Approximate number of characters in the SAM representation
    /// of this field
    int sam_length() const;

    /// Approximate number of bytes in the BAM representation of the
    /// SAM-formatted @a text
    /** Returns a conservative (i.e., generous) approximation of the
    BAM-formatted size of a (valid) SAM aux field of the given @a text
    of length @a text_length.
    FIXME Blah blah about returning 0 for invalid formatting.  */
    static int size_sam(const char* text, int text_length);

  private:
    // @cond private
    friend class alignment;
    friend char* format_sam(char* dest, const tagfield& aux);

    char tag_[2];
    char type_;
    char data[1];
    // @endcond
  };

  // @cond infrastructure
  class const_iterator;

  class iterator : public std::iterator<std::forward_iterator_tag, tagfield> {
  public:
    iterator() { }
    iterator(const iterator& it) : ptr(it.ptr) { }
    ~iterator() { }
    iterator& operator= (iterator it) { ptr = it.ptr; return *this; }

    bool operator== (iterator it) const { return ptr == it.ptr; }
    bool operator!= (iterator it) const { return ptr != it.ptr; }

    tagfield& operator* () const { return *reinterpret_cast<tagfield*>(ptr); }
    tagfield* operator-> () const { return reinterpret_cast<tagfield*>(ptr); }

    iterator& operator++ () { ptr += (*this)->size(); return *this; }
    iterator operator++ (int)
      { char* orig = ptr; ptr += (*this)->size(); return iterator(orig); }

  private:
    friend class alignment;
    friend class const_iterator;
    explicit iterator(char* p) : ptr(p) { }

    char* ptr;
  };

  class const_iterator :
    public std::iterator<std::forward_iterator_tag, tagfield,
			 std::ptrdiff_t, const tagfield*, const tagfield&> {
  public:
    const_iterator() { }
    const_iterator(const const_iterator& it) : ptr(it.ptr) { }
    const_iterator(alignment::iterator it) : ptr(it.ptr) { }
    ~const_iterator() { }
    const_iterator& operator= (const_iterator it)
      { ptr = it.ptr; return *this; }

    bool operator== (const_iterator it) const { return ptr == it.ptr; }
    bool operator!= (const_iterator it) const { return ptr != it.ptr; }

    const tagfield& operator* () const
      { return *reinterpret_cast<const tagfield*>(ptr); }

    const tagfield* operator-> () const
      { return reinterpret_cast<const tagfield*>(ptr); }

    const_iterator& operator++ () { ptr += (*this)->size(); return *this; }
    const_iterator operator++ (int)
      { const char* orig = ptr;
	ptr += (*this)->size(); return const_iterator(orig); }

  private:
    friend class alignment;
    friend std::ostream& operator<< (std::ostream&, const_iterator);
    explicit const_iterator(const char* p) : ptr(p) { }

    const char* ptr;
  };

  typedef iterator::reference reference;
  typedef const_iterator::reference const_reference;
  typedef iterator::difference_type difference_type;
  // @endcond

  iterator begin() { return iterator(p->auxen_data()); }
  const_iterator begin() const { return const_iterator(p->auxen_data()); }

  iterator end() { return iterator(p->end_data()); }
  const_iterator end() const { return const_iterator(p->end_data()); }

  iterator find(const char* tag);
  const_iterator find(const char* tag) const;

  bool empty() const { return p->auxen_data() == p->end_data(); }

  void push_back_sam(const char* aux_text, int aux_text_length);
  void push_back_sam(const char* aux_text)
    { push_back_sam(aux_text, strlen(aux_text)); }
  void push_back_sam(const std::string& aux_text)
    { push_back_sam(aux_text.c_str(), aux_text.length()); }

  template <typename ValueType>
  void push_back(const char* tag, ValueType value)
    { replace_(end(), end(), tag, value); }

  template <typename ValueType>
  iterator insert(iterator position, const char* tag, ValueType value)
    { return replace_(position, position, tag, value); }

  iterator erase(iterator position)
    { iterator next = position; return replace_gap(position, ++next, 0); }
  iterator erase(iterator start, iterator limit)
    { return replace_gap(start, limit, 0); }

  void clear() { replace_gap(begin(), end(), 0); }

  template <typename ValueType>
  iterator
  replace(iterator start, iterator limit, const char* tag, ValueType value)
    { return replace_(start, limit, tag, value); }
  //@}

  /** @name Field modifiers  */
  //@{
  void set_qname(const std::string& qname)
    { set_qname(qname.data(), qname.length()); }

  void set_flags(int flags);
  void set_rindex(int rindex);
  void set_rname(const std::string& rname);
  void set_pos(coord_t pos);
  void set_zpos(coord_t zpos);
  void set_mapq(int mapq);
  void set_cigar(const std::string& cigar) { set_cigar(cigar.c_str()); }
  void set_cigar(const std::vector<cigar_op>& cigar);
  void set_mate_rindex(int mate_rindex);
  void set_mate_rname(const std::string& mate_rname);
  void set_mate_pos(coord_t pos);
  void set_mate_zpos(coord_t zpos);
  void set_isize(scoord_t isize);
  void set_seq(const std::string& seq);

  /// Update an existing @a tag's value, or add a new auxiliary field
  template <typename ValueType>
  void set_aux(const char* tag, ValueType value)
    { iterator position = find(tag); iterator limit = position;
      if (position != end())  ++limit;
      replace_(position, limit, tag, value); }

  /// Update the existing auxiliary field's value
  template <typename ValueType>
  iterator set_aux(iterator position, ValueType value)
    {iterator next = position; return replace_(position, ++next, NULL, value);}

  /// Erase all auxiliary fields with the given @a tag
  /** @return The number of fields erased (usually no more than 1).  */
  int erase(const char* tag);
  //@}

  /** @name Additional field modifiers  */
  //@{
  void set_qname(const char* qname) { set_qname(qname, strlen(qname)); }

  void set_cigar(const char* cigar);

  // FIXME put length first?
  void set_raw_seq(const char* seq, int length);
  void set_raw_seq_qual(const char* seq, int length, const char* qual);
  //@}

  /// @name Derived information
  //@{
  /// The @e strand flag, as +1 or -1
  int strand() const { return (flags() & REVERSE_STRAND)? -1 : +1; }

  /// The @e strand flag, as '+' or '-'
  char strand_char() const { return (flags() & REVERSE_STRAND)? '-' : '+'; }

  /// The @e mate-strand flag, as +1 or -1
  int mate_strand() const { return (flags() & MATE_REVERSE_STRAND)? -1 : +1; }

  /// The @e mate-strand flag, as '+' or '-'
  /** Returns the mate read's strand information encoded in the flags field.
  This is only meaningful if the alignment in fact has a mate, i.e., if it
  is paired.
  @return Either '+' or '-'.  */
  char mate_strand_char() const
    { return (flags() & MATE_REVERSE_STRAND)? '-' : '+'; }

  /// The @e pair-order flags, as -1 or +1 (or 0 when unset)
  /** Returns the pairing order information encoded in the flags field.
  @return -1 when the @e first flag is set, +1 when the @e second flag is set,
  or 0 if neither (or both) is set.  */
  int order() const
    { return order_value[(flags() & (FIRST_IN_PAIR | SECOND_IN_PAIR)) >> 6]; }

  /// The number of reference bases spanned by the aligned read
  scoord_t cigar_span() const;
  // FIXME maybe name them cigar_qlength() and cigar_rlength() (was _span())?

  /// Rightmost position (1-based)
  coord_t right_pos() const  { return pos() + cigar_span() - 1; }
  /// Rightmost position (0-based)
  coord_t right_zpos() const { return zpos() + cigar_span() - 1; }
  //@}

  // FIXME nuke me or otherwise hide me!
  void dump_on(std::ostream&, const_iterator = const_iterator(0)) const;

  /// Pack sequence string into two-base-per-byte encoding
  /** Writes @a seq_length/2 (rounded up) bytes of encoded sequence data to
  @a dest.
  @param dest Output
  @param seq C-string
  @param seq_length blah */
  static void pack_seq(char* dest, const char* seq, int seq_length);

  // FIXME is this goofy?  maybe should just be a namespace-global function.
  /// Unpack two-base-per-byte encoded sequence data
  /** Assigns @a raw_seq's encoded sequence data to @a dest, as bases and
  potentially ambiguity characters, in uppercase.
  @param dest Blah output.
  @param raw_seq  Blah.  Not expected to be NUL-terminated.
  @param seq_length  Foo.  Note that this is the length of the read and the
  output string, @b not the length of the input @a raw_seq data, which is half
  as long.  */
  static void unpack_seq(std::string& dest, const char* raw_seq, int seq_length)
    { dest.resize(seq_length); unpack_seq(dest.begin(), raw_seq, seq_length); }

  /// Pack quality string into raw Phred encoding
  /** Writes @a seq_length bytes of encoded (by subtracting 33) quality data
  to @a dest.
  @param dest Output
  @param qual C-string
  @param seq_length blah */
  static void pack_qual(char* dest, const char* qual, int seq_length);

  /// Unpack raw-Phred-encoded quality string
  static void unpack_qual(std::string& dest, const char* phred, int seq_length)
    { dest.resize(seq_length); unpack_qual(dest.begin(), phred, seq_length); }

private:
  // @cond private
  friend class bamio;
  friend class samio;
  // FIXME Only friend because of cigar stuff and its unpack_seq/_qual usage
  friend char* format_sam(char*, const alignment&, const std::ios&);

  struct block_header {
    uint16_t capacity;
    uint16_t cindex;
  };

  struct bamcore {
    int32_t  rest_length;
    int32_t  rindex;
    int32_t  zpos;
    uint8_t  name_length;
    uint8_t  mapq;
    uint16_t bin;
    uint16_t cigar_length;
    uint16_t flags;
    int32_t  read_length;
    int32_t  mate_rindex;
    int32_t  mate_zpos;
    int32_t  isize;
  };

  class block {
  public:
    struct block_header h;
    struct bamcore c;
    char extra[1];

    // These measure the BAM-formatted payload; the actual capacity and
    // size of the alignment::block exceed these by sizeof(block_header).
    int capacity() const { return h.capacity; }
    int size() const { return sizeof(c.rest_length) + c.rest_length; }

    char* data() { return reinterpret_cast<char*>(&this->c); }
    char* end_data() { return data() + sizeof(c.rest_length) + c.rest_length; }

    char* name_data()  { return this->extra; }
    char* cigar_data() { return name_data() + c.name_length; }
    char* seq_data()   { return cigar_data() + sizeof(uint32_t)*c.cigar_length;}
    char* qual_data()  { return seq_data() + (c.read_length + 1) / 2; }
    char* auxen_data() { return qual_data() + c.read_length; }

    static block* create(int payload_size);
    static void destroy(block* block);
    static void copy(block* dest, const block* src);
  };

  block* p;

  void assign(int nfields, const std::vector<char*>& fields, int cindex);

  void sync() const { bin(); }

  void resize_unshare_copy(int payload_size);
  void resize_unshare_discard(int payload_size);

  void set_qname(const char* qname, int qname_length);

  const_iterator find_or_throw(const char* tag) const;

  char* replace_gap(char* start, char* limit, int gap_length);
  iterator replace_gap(iterator start, iterator limit, int gap_length)
    { return iterator(replace_gap(start.ptr, limit.ptr, gap_length)); }

  iterator replace_string(iterator start, iterator limit,
      const char* tag, char type, const char* value, int value_length);

  iterator replace_(iterator start, iterator limit,
		    const char* tag, const std::string& value)
    { return replace_string(start, limit, tag, 'Z',
			    value.data(), value.length()); }

  iterator replace_(iterator start, iterator limit,
		    const char* tag, const char* value)
    { return replace_string(start, limit, tag, 'Z', value, strlen(value)); }

  iterator replace_(iterator start, iterator limit, const char* tag,char value);
  iterator replace_(iterator start, iterator limit, const char* tag, int value);

  // TODO Add replace_ for float ('f'), double ('d'),
  // and const std::vector<uint8_t>& ('H') (and maybe a string adapter for 'H')

  iterator replace_(iterator start, iterator limit,
		    const char* tag, const_iterator value);

  // FIXME or are the char* ones public?
  static char* unpack_seq(char* dest, const char* raw_seq, int seq_length);
  static void unpack_seq(std::string::iterator dest,
			 const char* raw_seq, int seq_length);

  static char* unpack_qual(char* dest, const char* phred, int length);
  static void unpack_qual(std::string::iterator dest,
			  const char* phred, int length);

  static const uint16_t unknown_bin = 0xffff;
  static const int order_value[];
  static block empty_block;
  // @endcond
};

// @cond infrastructure
inline bool operator== (alignment::iterator it1, alignment::const_iterator it2)
  { return alignment::const_iterator(it1) == it2; }
inline bool operator!= (alignment::iterator it1, alignment::const_iterator it2)
  { return alignment::const_iterator(it1) != it2; }

template<> inline std::string alignment::cigar() const
  { std::string dest; return cigar(dest); }
template<> inline std::vector<cigar_op> alignment::cigar() const
  { std::vector<cigar_op> dest; return cigar(dest); }

template<> inline std::string alignment::tagfield::value() const
  { std::string dest; return value(dest); }

template<> const char* alignment::tagfield::value() const;
template<> int alignment::tagfield::value() const;
template<> char alignment::tagfield::value() const;
// @endcond

/// Compare alignments by genomic location
/** The alignments are ordered by reference index (with unmapped, -1, sorting
last), position, query name, and ordering flag.
@relatesalso alignment */
bool operator< (const alignment& a, const alignment& b);

/// Compare alignments by genomic location, similarly to @c operator< above
/** @relatesalso alignment */
inline bool operator> (const alignment& a, const alignment& b) { return b < a; }

/// Compare alignments by query name
/** @return Less than, equal to, or greater than 0, similarly to @c strcmp().
@relatesalso alignment */
inline int cmp_by_qname(const alignment& a, const alignment& b)
  { return strcmp(a.qname_c_str(), b.qname_c_str()); }

/// Swap two alignments
/** @relatesalso alignment */
inline void swap(alignment& a, alignment& b) { a.swap(b); }

#if 0
/** @brief Read an alignment from the stream
@details Extracts a single SAM-formatted alignment record from the input
@a stream into @a aln.  Sets @c failbit if the first line available is not
an alignment record, for example if it starts with an '@' character.
@relatesalso alignment */
// FIXME NUKE-ME  What on earth do we do for a sam::collection?
// Could hook it up to one of the istream's pword()s, but that's insanely OTT.
std::istream& operator>> (std::istream& stream, alignment& aln);
#endif

/// Print an alignment to the stream
/** Writes an alignment record to a stream as text in SAM format, @b without
a trailing newline character.
@relatesalso alignment */
std::ostream& operator<< (std::ostream& stream, const alignment& aln);

/// Write an alignment to @a dest in SAM format
/** @param dest  Character array to be written to; must have space for at least
alignment::sam_length() characters.
@param aln  The alignment record to be formatted.
@param format  Format flags controlling how the record (in particular,
its @em FLAG field) should be formatted.
@return  A pointer to the first unused character position in @a dest.
@relatesalso alignment */
char* format_sam(char* dest, const alignment& aln, const std::ios& format);

/// Write alignment flags to @a dest in SAM format
/** @param dest  Character array to be written to; must have space for
at least 16 characters.
@param flags  The alignment flags to be formatted.
@param format  The representation written will be symbolic if @a format.flags()
has @c boolalpha set; otherwise numeric as determined by @c basefield,
but always in uppercase irrespective of @c uppercase.
@return  A pointer to the first unused character position in @a dest.
@relatesalso alignment */
char* format_sam(char* dest, int flags, const std::ios& format);

/// Print an auxiliary field to the stream
/** Writes an auxiliary field to a stream as text in SAM format.
@relatesalso alignment::tagfield */
std::ostream& operator<< (std::ostream& stream, const alignment::tagfield& aux);

/// Write an auxiliary field to @a dest in SAM format
/** @return  A pointer to the first unused character position in @a dest.
@relatesalso alignment::tagfield */
char* format_sam(char* dest, const alignment::tagfield& aux);

/// Print an iterator or const_iterator to the stream
std::ostream& operator<< (std::ostream& stream, alignment::const_iterator it);

} // namespace sam

#endif
