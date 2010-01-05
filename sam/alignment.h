/** @file   sam/alignment.h
    @brief  Classes and functions for SAM/BAM alignment records
*/

#ifndef CANSAM_ALIGNMENT_H
#define CANSAM_ALIGNMENT_H

#include <string>
#include <vector>
//#include <algorithm>  // FIXME NUKE-ME, but for specialising std::swap IIRC
#include <iterator>

#include <ostream> // FIXME NUKE-ME, but figure out who's getting us <iosfwd>

#include <stdint.h>

#include "sam/types.h"
#include "sam/collection.h"

namespace sam {

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
  DUPLICATE           = 0x400  ///< PCR duplicate or optical duplicate
};

/** Returns the BAM bin number for an alignment spanning [@a pos, @a right],
    i.e., a 1-based range.  */
int calc_bin(coord_t pos, coord_t right);

/** Returns the BAM bin number for an alignment spanning [@a zpos, @a zright],
    i.e., a 0-based range.  */
int calc_zbin(coord_t zpos, coord_t zright);

/** @class sam::alignment sam/alignment.h
    @brief SAM/BAM alignment record */
class alignment {
public:
  /// Construct an empty alignment.
  alignment() : p(&empty) { }

#if 0
  /// Construct an alignment by splitting up a tab-separated text string.
  // FIXME what reference_thingie does this use?
  explicit alignment(const std::string& line) { assign(line); }
#endif

  /// Construct a copy of an alignment.
  alignment(const alignment& aln);

  /// Destroy this alignment object.
  ~alignment() { if (p != &empty)  block::destroy(p); }

  /// Copy an alignment.
  alignment& operator= (const alignment& aln);

  /// Assign to this alignment by splitting up a tab-separated text string.
  // FIXME what reference_thingie does this use?
  alignment& assign(const std::string& line);

  /// Swap this alignment with another.
  void swap(alignment& aln) { block* tmp = p; p = aln.p; aln.p = tmp; }

  void reserve(int auxsize);

  /** @name Field accessors */
  //@{
  /// Query name
  std::string qname() const
    { return std::string(p->data() + p->name_offset(), p->c.name_length - 1); }

  /// Bitwise combination of alignment_flags
  int flags() const { return p->c.flags; }

  int rindex() const { return p->c.rindex; } ///< Reference identifier (or -1)

  /// Reference name (or '*' @c *)
  std::string rname() const
    ;// FIXME { return collection::find(p->h.cindex).rname(p->c.rindex); }

  coord_t pos() const  { return p->c.zpos+1; } ///< Leftmost position (1-based)
  coord_t zpos() const { return p->c.zpos; }   ///< Leftmost position (0-based)

  /// BAM bin number
  int bin() const
    { if (p->c.bin == unknown_bin)  p->c.bin = calc_zbin(zpos(), right_zpos());
      return p->c.bin; }

  int mapq() const { return p->c.mapq; } ///< Mapping quality

  std::string cigar() const; ///< Cigar string
  // FIXME raw cigar

  /// Mate's reference identifier (or -1; for paired reads)
  int mate_rindex() const { return p->c.mate_rindex; }
  std::string mate_rname() const; ///< Mate's reference name (for paired reads)

  /// Mate's leftmost position (1-based)
  coord_t mate_pos() const  { return p->c.mate_zpos+1; }
  /// Mate's leftmost position (0-based)
  coord_t mate_zpos() const { return p->c.mate_zpos; }

  scoord_t isize() const { return p->c.isize; } ///< Insert size

  int length() const { return p->c.read_length; } ///< Read length

  /// Sequence
  std::string seq() const
    { std::string dest(length(), 'X');
      unpack_seq(dest.begin(), raw_seq(), length());
      return dest; }

  /// Sequence (packed as two bases per byte)
  const char* raw_seq() const { return p->data() + p->seq_offset(); }

  /// Quality string
  // FIXME uh, no, needs offsetted
  std::string qual() const
    { return std::string(p->data() + p->qual_offset(), p->c.read_length); }

  /// Value of the auxiliary field with the given @a tag
  std::string aux(const char* tag) const;
  std::string aux(const char* tag, const std::string& default_value) const;

  int aux_int(const char* tag) const;
  int aux_int(const char* tag, int default_value) const;

  float aux_float(const char* tag) const;
  float aux_float(const char* tag, float default_value) const;

  double aux_double(const char* tag) const;
  double aux_double(const char* tag, double default_value) const;
  //@}

  /** @name Auxiliary fields as a collection
  Alignment records provide limited collection-style access to their
  auxiliary fields.

  The @c iterator and @c const_iterator classes are <b>forward iterators</b>
  providing all the usual iterator functionality: copying, assignment,
  pre- and post-increment, equality and inequality tests, and dereferencing,
  which produces an aux_field through which the pointed-to auxiliary field's
  properties can be accessed.

  This access is read-only; to change a field's value, use the iterator
  variants of set_aux() listed below.  */
  //@{
  /** @class sam::alignment::aux_field sam/alignment.h
      @brief Helper class representing an auxiliary field as seen via an
	     iterator

  Dereferencing a sam::alignment @c iterator or @c const_iterator produces
  (a reference to) an instance of this class.

  @note There are no mutator methods, even if it is a mutable @c iterator that
  has been dereferenced; use the iterator variants of sam::alignment::set_aux()
  to change the value of an auxiliary field via an @c iterator.  */
  // FIXME maybe rename to tagfield, and there will also be a header::tagfield
  class aux_field {
  public:
    /// Two-character field tag
    std::string tag() const { return std::string(tag_, sizeof tag_); }

    /// Field type, as it would appear in a BAM file
    /** @return One of 'A', 'Z', 'i', 'c', 'S', etc, including the subtypes
    that only appear in a BAM file.  */
    char type() const { return type_; }

    /// Field value, as it would appear in a SAM file
    /** @return The text string representation of the field's value (rather
    than any binary representation that might appear in a BAM file).  */
    std::string value() const;

    /// Field value as an integer
    /// (or throws if this field's type is non-integral)
    int value_int() const;

    // TODO  Implement value_coord(), value_float(), value_double()

    /// Returns whether the field's tag is the same as the given @a key_tag
    bool tag_equals(const char* key_tag) const
      { return tag_[0] == key_tag[0] && tag_[1] == key_tag[1]; }

    /// Number of bytes in the BAM representation of this field
    int size() const;

  private:
    char tag_[2];
    char type_;
    char data[1];
  };

  // @cond infrastructure
  class const_iterator;

  class iterator : public std::iterator<std::forward_iterator_tag, aux_field> {
  public:
    iterator() { }
    iterator(const iterator& it) : aux(it.aux) { }
    ~iterator() { }
    iterator& operator= (iterator it) { aux = it.aux; return *this; }

    aux_field& operator* () const { return *reinterpret_cast<aux_field*>(aux); }
    aux_field* operator-> () const { return reinterpret_cast<aux_field*>(aux); }

    iterator& operator++ () { aux += (*this)->size(); return *this; }
    iterator operator++ (int)
      { char* orig = aux; aux += (*this)->size(); return iterator(orig); }

  private:
    friend class alignment;
    friend class const_iterator;
    explicit iterator(char* p) : aux(p) { }

    char* aux;
  };

  class const_iterator :
    public std::iterator<std::forward_iterator_tag, aux_field,
			 ptrdiff_t, const aux_field*, const aux_field&> {
  public:
    const_iterator() { }
    const_iterator(const const_iterator& it) : aux(it.aux) { }
    const_iterator(iterator it) : aux(it.aux) { }
    ~const_iterator() { }
    const_iterator& operator= (const_iterator it)
      { aux = it.aux; return *this; }

    const aux_field& operator* () const
      { return *reinterpret_cast<const aux_field*>(aux); }

    const aux_field* operator-> () const
      { return reinterpret_cast<const aux_field*>(aux); }

    const_iterator& operator++ () { aux += (*this)->size(); return *this; }
    const_iterator operator++ (int)
      { const char* orig = aux;
	aux += (*this)->size(); return const_iterator(orig); }

    friend inline bool operator== (const_iterator it1, const_iterator it2)
      { return it1.aux == it2.aux; }
    friend inline bool operator!= (const_iterator it1, const_iterator it2)
      { return it1.aux != it2.aux; }

    friend std::ostream& operator<< (std::ostream& stream, const_iterator it);

  private:
    friend class alignment;
    explicit const_iterator(const char* p) : aux(p) { }

    const char* aux;
  };

  typedef iterator::reference reference;
  typedef const_iterator::reference const_reference;
  typedef iterator::difference_type difference_type;
  // @endcond

  iterator begin() { return iterator(p->data() + p->auxen_offset()); }
  const_iterator begin() const
    { return const_iterator(p->data() + p->auxen_offset()); }

  iterator end() { return iterator(p->data() + p->end_offset()); }
  const_iterator end() const
    { return const_iterator(p->data() + p->end_offset()); }

  iterator find(const char* tag);
  const_iterator find(const char* tag) const;

  void push_back(const char* tag, const std::string& value)
    { replace(end(), end(), tag, value); }
  void push_back(const char* tag, int value)
    { replace(end(), end(), tag, value); }

#if 1
  iterator insert(iterator position, const char* tag, const std::string& value)
    { return replace(position, position, tag, value); }
  iterator insert(iterator position, const char* tag, int value)
    { return replace(position, position, tag, value); }
#else
  template<typename value_type>
  iterator insert(iterator position, const char* tag, value_type value)
    { return replace(position, position, tag, value); }
#endif

  iterator erase(iterator position)
    //{ return replace_gap(position, next(position), 0); }
    { iterator first = position; return replace_gap(first, ++position, 0); }
  iterator erase(iterator start, iterator limit)
    { return replace_gap(start, limit, 0); }

  void clear() { replace_gap(begin(), end(), 0); }

  iterator replace(iterator start, iterator limit,
		   const char* tag, const std::string& value);
  iterator replace(iterator start, iterator limit, const char* tag, int value);
  //@}

  /** @name Additional field accessors
      Accessors returning strings also have corresponding versions returning
      C strings, which are slightly cheaper as generally a NUL-terminated
      string is already available and no string construction memory allocation
      is needed.
      
      Returned pointers point within the blahblah, and are valid as long as
      no non-const member functions (including the destructor!) are called
      for this alignment.  FIXME  */
  //@{
  /// Query name
  const char* qname_c_str() const { return p->data() + p->name_offset(); }

  /// Assigns query name to @a dest (and returns @a dest)
  std::string& qname(std::string& dest) const
    { return dest.assign(p->data() + p->name_offset(), p->c.name_length - 1); }

  const char* rname_c_str() const;      ///< Reference name (or @c NULL)
  const char* mate_rname_c_str() const; ///< Mate's reference name (or @c NULL)

  /// Assigns sequence to @a dest (and returns @a dest)
  std::string& seq(std::string& dest) const
    { unpack_seq(dest, raw_seq(), length()); return dest; }

  /// Quality string (@b not NUL-terminated)
  const char* qual_data() const { return p->data() + p->qual_offset(); }

  const char* aux_c_str(const char* tag) const;
  //@}

  /** @name Field modifiers  */
  //@{
  void set_qname(const std::string& qname);
  void set_flags(int flags);
  void set_rindex(int rindex);
  void set_rname(const std::string& rname);
  void set_pos(coord_t pos);
  void set_zpos(coord_t zpos);
  void set_mapq(int mapq);
  void set_cigar(const std::string& cigar);
  void set_mate_rindex(int mate_rindex);
  void set_mate_rname(const std::string& mate_rname);
  void set_mate_pos(coord_t pos);
  void set_mate_zpos(coord_t zpos);
  void set_isize(scoord_t isize);

  // FIXME put length first?
  void set_seq(const std::string& seq);
  void set_raw_seq(const char* seq, int length);
  void set_raw_seq_qual(const char* seq, int length, const char* qual);

  /* Replace() primitives:
  replace_gap(iterator start, iterator limit, int gap_length); [returns rewired start]
  replace(iterator start, iterator limit, const char* tag, const std::string& value);
    { replace_gap(st,lim,3+value.length()+1); memcpy; }
  replace(iterator start, iterator limit, const char* tag, int value);
    { if (value <= INT8_MAX) { replace_gap(st,lim,3+1); cvt_into; } else if... }
  // etc

  Now implement in terms of this:

  void set_aux(const char* tag, const std::string& value);	// lim = pos = find(tag); if (pos != end()) ++lim; replace(pos, lim, tag, value)
  void set_aux(iterator position, const std::string& value);	// replace(pos, next(pos), orig->tag(), value)

  void push_back(const char* tag, const std::string& value);	// replace(end(), end(), tag, value)
  iterator insert(iterator position, const char* tag, value);	// replace(pos, pos, tag, value)

  void erase(const char* tag);					// pos = find(tag); if (pos != end()) replace_gap(pos, next(pos), 0)
  iterator erase(iterator position)				// replace_gap(pos, next(pos), 0)
  iterator erase(iterator first, iterator last);		// replace_gap(first, last, 0)

  void clear() { erase(begin(), end()); }

  */

  /// Update an existing @a tag's value, or add a new auxiliary field
  void set_aux(const char* tag, const std::string& value);
  void set_aux(const char* tag, int value);

  iterator set_aux(iterator position, const std::string& value);
  iterator set_aux(iterator position, int value);

  void erase(const char* tag);
  //@}

  /// @name Derived information
  //@{
  /// The @e strand flag, as '+' or '-'
  /** Returns the strand information encoded in the flags field.
  @return Either '+' or '-'.  */
  char strand() const { return (flags() & REVERSE_STRAND)? '-' : '+'; }

  /// @brief The @e mate-strand flag, as '+' or '-'.
  /// @details Returns the mate read's strand information encoded in the
  /// flags field.  This is only meaningful if the alignment in fact has
  /// a mate, i.e., if it is paired.
  /// @return Either '+' or '-'.
  char mate_strand() const {return (flags() & MATE_REVERSE_STRAND)? '-' : '+';}

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

  // FIXME maybe should be private/friend?
  int approx_sam_record_length() const;
  void sam_record(char*, int) const;
  //@}

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

  // FIXME probably will end up as  std::string& dest
  static void unpack_qual(char* dest, const char* phred, int seq_length);

private:
  // @cond private
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

    int capacity() const { return h.capacity; }
    int size() const
      { return sizeof(block_header) + sizeof(c.rest_length) + c.rest_length; }

    char* data() { return reinterpret_cast<char*>(&this->c); }

    int name_offset() const  { return sizeof(bamcore); }
    int cigar_offset() const { return name_offset() + c.name_length; }
    int seq_offset() const   { return cigar_offset() + 4 * c.cigar_length; }
    int qual_offset() const  { return seq_offset() + (c.read_length + 1) / 2; }
    int auxen_offset() const { return qual_offset() + c.read_length; }
    int end_offset() const   { return sizeof(c.rest_length) + c.rest_length; }

    static block* create(int capacity);
    static void destroy(block* block);
    static void copy(block* dest, const block* src);
  };

  block* p;

public: // FIXME figure out how to get samstream_base::samio access to this
  void assign(int nfields, const std::vector<char*>& fields /*, refthingie*/);
private:

  void sync() const;
  //FIXME NUKEME void destroy();  // FIXME prob should be const (!)

  void resize_unshare_copy(int size);
  void resize_unshare_discard(int size);

  iterator replace_gap(iterator start, iterator limit, int gap_length);

  static void unpack_seq(std::string::iterator dest,
			 const char* raw_seq, int seq_length);

  static const uint16_t unknown_bin = 0xffff;
  static const int order_value[];
  static block empty;
  // @endcond
};

/// Compare alignments by genomic location
/** The alignments are ordered by reference index (with unmapped, -1, sorting
last), position, query name, and ordering flag.
@relatesalso alignment */
bool operator< (const alignment& a, const alignment& b);

/// Compare alignments by genomic location, similarly to @c operator< above
/** @relatesalso alignment */
inline bool operator> (const alignment& a, const alignment& b) { return b < a; }

/// Swap two alignments
/** @relatesalso alignment */
inline void swap(alignment& a, alignment& b) { a.swap(b); }

/** @brief Read an alignment from the stream
@details Extracts a single SAM-formatted alignment record from the input
@a stream into @a aln.  Sets @c failbit if the first line available is not
an alignment record, for example if it starts with an '@' character.
@relatesalso alignment */
// FIXME NUKE-ME  What on earth do we do for a sam::collection?
// Could hook it up to one of the istream's pword()s, but that's insanely OTT.
std::istream& operator>> (std::istream& stream, alignment& aln);

/// Print an alignment to the stream
/** Writes an alignment record to a %stream as text in SAM format, @b without
a trailing newline character.
@relatesalso alignment */
std::ostream& operator<< (std::ostream& stream, const alignment& aln);

} // namespace sam

#endif
