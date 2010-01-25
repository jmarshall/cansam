/** @file   sam/alignment.h
    @brief  Classes and functions for SAM/BAM alignment records
*/

#ifndef CANSAM_ALIGNMENT_H
#define CANSAM_ALIGNMENT_H

#include <string>
#include <vector>
#include <iterator>
//#include <algorithm>  // FIXME NUKE-ME, but for specialising std::swap IIRC
#include <iosfwd>

#include <stdint.h>

#include "sam/types.h"
#include "sam/header.h"
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

/// Returns the BAM bin number (1-based)
/** Returns the BAM bin number for an alignment spanning [@a pos, @a right],
    i.e., a 1-based range.  */
int calc_bin(coord_t pos, coord_t right);

/// Returns the BAM bin number (0-based)
/** Returns the BAM bin number for an alignment spanning [@a zpos, @a zright],
    i.e., a 0-based range.  */
int calc_zbin(coord_t zpos, coord_t zright);

/** @class sam::alignment sam/alignment.h
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
  - @c const_iterator
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

  void reserve(int auxsize); // FIXME How much should @a size measure?

  /** @name Field accessors
  Two variants are provided for the @em POS and @em MPOS fields: @c %pos()
  et al return 1-based coordinates, while @c %zpos() et al return the same
  positions in 0-based coordinates.  */
  //@{
  std::string qname() const
    { return std::string(p->data() + p->name_offset(), p->c.name_length - 1); }

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
  std::string cigar() const;
  // TODO raw cigar

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

  // FIXME uh, no, needs offsetted
  std::string qual() const
    { return std::string(p->data() + p->qual_offset(), p->c.read_length); }

  int length() const { return p->c.read_length; } ///< Read length

  /// BAM bin number (derived, if necessary, from @em POS and @em CIGAR)
  int bin() const
    { if (p->c.bin == unknown_bin)  p->c.bin = calc_zbin(zpos(), right_zpos());
      return p->c.bin; }

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

  /** @name Additional field accessors
  Some accessors returning strings also have corresponding versions returning
  C strings or (non-NUL-terminated) @c char arrays.  These versions are cheaper
  as a suitable representation already exists within the alignment object, so
  no @c std::string construction or memory allocation is necessary.

  Pointers returned by an alignment object become invalid when any of
  that object's non-const member functions are subsequently called.  */
  //@{
  /// Query name
  const char* qname_c_str() const { return p->data() + p->name_offset(); }

  /// Assigns query name to @a dest (and returns @a dest)
  std::string& qname(std::string& dest) const
    { return dest.assign(p->data() + p->name_offset(), p->c.name_length - 1); }

  ///< Reference name (or @c NULL)
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
  const char* seq_raw_data() const { return p->data() + p->seq_offset(); }

  /// Quality string (not NUL-terminated)
  const char* qual_data() const { return p->data() + p->qual_offset(); }

  const char* aux_c_str(const char* tag) const;
  //@}

  /** @name Auxiliary fields as a collection
  Alignment records provide limited collection-style access to their
  auxiliary fields.

  The @c sam::alignment::iterator and @c sam::alignment::const_iterator classes
  are <b>forward iterators</b> providing all the usual iterator functionality:
  copying, assignment, pre- and post-increment, equality and inequality tests,
  and dereferencing, which produces an alignment::tagfield through which the
  pointed-to auxiliary field's properties can be accessed (but not modified;
  see also the iterator variant of set_aux() below).

  The possible types for the @c ValueType parameter are listed below.  */
  //@{
  /** @class sam::alignment::tagfield sam/alignment.h
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
    friend class alignment;

    char tag_[2];
    char type_;
    char data[1];
  };

  // @cond infrastructure
  typedef std::char_traits<char> traits_type;
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
			 ptrdiff_t, const tagfield*, const tagfield&> {
  public:
    const_iterator() { }
    const_iterator(const const_iterator& it) : ptr(it.ptr) { }
    const_iterator(iterator it) : ptr(it.ptr) { }
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

  iterator begin() { return iterator(p->data() + p->auxen_offset()); }
  const_iterator begin() const
    { return const_iterator(p->data() + p->auxen_offset()); }

  iterator end() { return iterator(p->data() + p->end_offset()); }
  const_iterator end() const
    { return const_iterator(p->data() + p->end_offset()); }

  iterator find(const char* tag);
  const_iterator find(const char* tag) const;

  bool empty() const { return p->auxen_offset() == p->end_offset(); }

  template <typename ValueType>
  void push_back(const char* tag, ValueType value)
    { replace_(end(), end(), tag, value); }

  template <typename ValueType>
  iterator insert(iterator position, const char* tag, ValueType value)
    { return replace_(position, position, tag, value); }

  iterator erase(iterator position)
    { return replace_gap(position, next(position), 0); }
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
  void set_cigar(const std::string& cigar);
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
    { return replace_(position, next(position), NULL, value); }

  /// Erase all auxiliary fields with the given @a tag
  /** @return The number of fields erased (usually no more than 1).  */
  int erase(const char* tag);
  //@}

  /** @name Additional field modifiers  */
  //@{
  void set_qname(const char* qname)
    { set_qname(qname, traits_type::length(qname)); }

  // FIXME put length first?
  void set_raw_seq(const char* seq, int length);
  void set_raw_seq_qual(const char* seq, int length, const char* qual);
  //@}

  /// @name Derived information
  //@{
  /// The @e strand flag, as '+' or '-'
  /** Returns the strand information encoded in the flags field.
  @return Either '+' or '-'.  */
  char strand() const { return (flags() & REVERSE_STRAND)? '-' : '+'; }

  /// The @e mate-strand flag, as '+' or '-'
  /** Returns the mate read's strand information encoded in the flags field.
  This is only meaningful if the alignment in fact has a mate, i.e., if it
  is paired.
  @return Either '+' or '-'.  */
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

  void set_qname(const char* qname, int qname_length);

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
    { return replace_string(start, limit, tag, 'Z',
			    value, traits_type::length(value)); }

  iterator replace_(iterator start, iterator limit, const char* tag, int value);

  // TODO Add replace_ for char ('A'), float ('f'), double ('d'),
  // and const std::vector<uint8_t>& ('H') (and maybe a string adapter for 'H')

  iterator replace_(iterator start, iterator limit,
		    const char* tag, const_iterator value);

  static void unpack_seq(std::string::iterator dest,
			 const char* raw_seq, int seq_length);

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
// @endcond

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
/** Writes an alignment record to a stream as text in SAM format, @b without
a trailing newline character.
@relatesalso alignment */
std::ostream& operator<< (std::ostream& stream, const alignment& aln);

/// Print an auxiliary field to the stream
/** Writes an auxiliary field to a stream as text in SAM format.
@relatesalso alignment::tagfield */
std::ostream& operator<< (std::ostream& stream, const alignment::tagfield& aux);

/// Print an iterator or const_iterator to the stream
std::ostream& operator<< (std::ostream& stream, alignment::const_iterator it);

} // namespace sam

#endif
