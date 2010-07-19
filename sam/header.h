/** @file   sam/header.h
    @brief  Classes for SAM/BAM headers
*/

#ifndef CANSAM_HEADER_H
#define CANSAM_HEADER_H

#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <iosfwd>

#include "sam/types.h"

namespace sam {

class sambamio;
class bamio;
class samio;
class collection;

/** @class sam::header sam/header.h
    @brief SAM/BAM header record, representing a single '@@' header line

This class blah FIXME blah container-style blah, so a header object itself
is amenable to the usual container operations, such as iterating through the
fields with @c begin() and @c end(), and modifying them with @c erase(),
@c push_back(), etc.

A SAM header probably ought not to have two fields with the same tag,
though this is not explicitly stated in the specification and is not enforced
by this class.  Thus find() and field() in fact return the @b first field with
the specified tag.  In the unlikely event that you expect a tag to occur in
several fields, you should iterate until @c end() checking for other instances
of the tag yourself.  */
class header {
public:
  header() { }
  // FIXME or should be protected
  // FIXME Hmmm... who ate all the tabs?  or ought to?
  explicit header(const std::string& line) : str_(line), cstr_(str_.c_str()) { }
  virtual ~header() { }

  // FIXME Infrastructure: ctors, copy ctor, assign, swap, reserve
  // FIXME Maybe... or maybe not sensible for polymorphic heap class

  /// The header's two-character type code
  std::string type() const;

  /// Returns whether this header's type code is the same as the
  /// given @a key_type
  bool type_equals(const char* key_type) const
    { return cstr_[0]=='@' && cstr_[1]==key_type[0] && cstr_[2]==key_type[1]; }

  /// Returns the (tab-separated) string representing this header
  std::string str() const;

  /// Returns the (exact) number of characters in the SAM representation
  /// of this header
  int sam_length() const { return str_.length(); }

  /// Returns the @a tag's value, or throws an exception if not found
  /** Returns the value of the header field with the given @a tag,
  or throws sam::exception if there is no such field or sam::bad_format
  if the field cannot be expressed as a @a ValueType.
  @note In most cases, it will be necessary to use explicit syntax such
  as <tt>field<int>("LN")</tt> to specify which type is desired.  */
  template <typename ValueType>
  ValueType field(const char* tag) const
    { return find_or_throw(tag)->template value<ValueType>(); }

  /// Returns the @a tag's value, or @a default_value if not found
  /** Returns the value of the header field with the given @a tag,
  or @a default_value if there is no such field; or throws sam::bad_format
  if the field cannot be expressed as a @a ValueType.  */
  template <typename ValueType>
  ValueType field(const char* tag, ValueType default_value) const
    { const_iterator it = find(tag);
      return (it != end())? it->value<ValueType>() : default_value; }

  /** @name Additional field accessors
  These field accessors operate similarly to those above; but, rather than
  returning a newly-constructed string, they assign to a destination string
  and return that string.  */
  //@{
  std::string& field(std::string& dest, const char* tag) const
    { return find_or_throw(tag)->value(dest); }

  std::string& field(std::string& dest,
		     const char* tag, const char* default_value) const
    { const_iterator it = find(tag);
      return (it != end())? it->value(dest) : dest.assign(default_value); }

  std::string& field(std::string& dest,
		     const char* tag, const std::string& default_value) const
    { const_iterator it = find(tag);
      return (it != end())? it->value(dest) : dest.assign(default_value); }
  //@}

  /** @name Container functionality
  Headers that have fields provide limited container-style access
  to their fields.

  The @c sam::header::iterator and @c sam::header::const_iterator classes are
  <b>bidirectional iterators</b> providing all the usual iterator functionality:
  copying, assignment, pre- and post-increment and -decrement, equality and
  inequality tests, and dereferencing, which produces a header::tagfield
  through which the pointed-to field's properties can be accessed (but not
  modified; see also the iterator variant of set_field() below).  */
  //@{
  class const_iterator;

  /** @class sam::header::tagfield sam/header.h
      @brief Helper class representing a header field as seen via an iterator

  Dereferencing a sam::header @c iterator or @c const_iterator produces
  (a reference to) an instance of this class.

  @note There are no mutator methods, even if it is a mutable @c iterator that
  has been dereferenced; use sam::header::set_field() to change the value of
  a header field via an @c iterator.  */
  class tagfield {
  public:
    /// Two-character field tag
    std::string tag() const;

    /// Field value
    template <typename ValueType>
    ValueType value() const;

    /// Assigns field value to @a dest (and returns @a dest)
    std::string& value(std::string& dest) const { return dest.assign(data_); }

    /// Returns whether the field's tag is the same as the given @a key_tag
    bool tag_equals(const char* key_tag) const
      {return tag_[0] == key_tag[0] && tag_[1] == key_tag[1] && colon_ == ':';}

  private:
    friend class const_iterator;
    friend std::ostream& operator<< (std::ostream&, const tagfield&);

    inline static const char* next(const char* s)
      { while (*++s != '\0') { }  return s; }

    inline static const char* prev(const char* s)
      { while (*--s != '\0') { }  return s; }

    char nul_;
    char tag_[2];
    char colon_;
    char data_[1];
  };

  // @cond infrastructure
  typedef std::char_traits<char> traits_type;

  class const_iterator :
    public std::iterator<std::bidirectional_iterator_tag, tagfield,
			 ptrdiff_t, const tagfield*, const tagfield&> {
  public:
    const_iterator() { }
    const_iterator(const const_iterator& it) : ptr(it.ptr) { }
    ~const_iterator() { }
    const_iterator& operator= (const_iterator it) {ptr = it.ptr; return *this;}

    bool operator== (const_iterator it) const { return ptr == it.ptr; }
    bool operator!= (const_iterator it) const { return ptr != it.ptr; }

    const tagfield* operator-> () const
      { return reinterpret_cast<const tagfield*>(ptr); }
    const tagfield& operator* () const { return *operator->(); }

    const_iterator& operator++ () { ptr = tagfield::next(ptr); return *this; }
    const_iterator operator++ (int)
      { const char* orig = ptr; ++(*this); return const_iterator(orig); }

    const_iterator& operator-- () { ptr = tagfield::prev(ptr); return *this; }
    const_iterator operator-- (int)
      { const char* orig = ptr; --(*this); return const_iterator(orig); }

  protected:
    friend class header;
    friend std::ostream& operator<< (std::ostream&, const_iterator);
    explicit const_iterator(const char* p) : ptr(p) { }

    const char* ptr;
  };

  class iterator : public const_iterator {
  public:
    iterator() { }
    iterator(const iterator& it) : const_iterator(it.ptr) { }
    ~iterator() { }
    iterator& operator= (iterator it) { ptr = it.ptr; return *this; }

    iterator& operator++ () { const_iterator::operator++(); return *this; }
    iterator operator++ (int)
      { const char* orig = ptr; ++(*this); return iterator(orig); }

    iterator& operator-- () { const_iterator::operator--(); return *this; }
    iterator operator-- (int)
      { const char* orig = ptr; --(*this); return iterator(orig); }

  private:
    friend class header;
    explicit iterator(const char* p) : const_iterator(p) { }
  };

  typedef iterator::reference reference;
  typedef const_iterator::reference const_reference;
  typedef const_iterator::difference_type difference_type;
  // @endcond

  iterator begin() { return iterator(cstr_ + 3); }
  const_iterator begin() const { return const_iterator(cstr_ + 3); }

  iterator end() { return iterator(cstr_ + str_.length()); }
  const_iterator end() const { return const_iterator(cstr_ + str_.length()); }

  iterator find(const char* tag) { return iterator(cstr_ + find_or_eos(tag)); }
  const_iterator find(const char* tag) const
    { return const_iterator(cstr_ + find_or_eos(tag)); }

  bool empty() const { return str_.length() <= 3; }

  template <typename ValueType>
  void push_back(const char* tag, ValueType value)
    { replace_(str_.length(), 0, tag, value); }

  template <typename ValueType>
  iterator insert(iterator position, const char* tag, ValueType value)
    { return replace_(position.ptr - cstr_, 0, tag, value); }

  iterator erase(iterator position)
    { iterator next = position; return erase(position, ++next); }

  iterator erase(iterator start, iterator limit);
  void clear();

  template <typename ValueType>
  iterator
  replace(iterator start, iterator limit, const char* tag, ValueType value)
    { return replace_(start.ptr - cstr_, limit.ptr - start.ptr, tag, value); }
  //@}

  /// Update an existing @a tag's value, or add a new header field
  template <typename ValueType>
  void set_field(const char* tag, ValueType value);

  /// Update the existing header field's value
  template <typename ValueType>
  iterator set_field(iterator position, ValueType value)
    { iterator next = position; return replace(position, ++next, NULL, value); }

  /// Erase all fields with the given @a tag
  /** @return The number of fields erased (usually no more than 1).  */
  int erase(const char* tag);

protected:
  /// Enable derived classes to update their state when a header is modified
  /** Called whenever a header is modified.  Blah blah blah.
  @note Derived classes augmenting this method should ensure that the first
  thing their overriding function does is invoke @c header::sync().  */
  virtual void sync() { cstr_ = str_.c_str(); }

private:
  // @cond private
  friend char* format_sam(char* dest, const header& header);
  // @endcond

  size_t find_or_eos(const char* tag) const;
  const_iterator find_or_throw(const char* tag) const;

  iterator replace_string(size_t pos, size_t length,
			  const char* tag, const char* value, int value_length);

  iterator replace_(size_t pos, size_t length,
		    const char* tag, const std::string& value)
    { return replace_string(pos, length, tag, value.data(), value.length()); }

  iterator replace_(size_t pos, size_t length,
		    const char* tag, const char* value)
    {return replace_string(pos, length, tag, value, traits_type::length(value));}

  iterator replace_(size_t pos, size_t length, const char* tag, int value);
  iterator replace_(size_t pos, size_t length,
		    const char* tag, const_iterator value);

  std::string str_;
  const char* cstr_;
};

template<> inline std::string header::tagfield::value() const
  { return std::string(data_); }
template<> inline const char* header::tagfield::value() const
  { return data_; }

template<> int header::tagfield::value() const;
template<> coord_t header::tagfield::value() const;

template <typename ValueType>
void header::set_field(const char* tag, ValueType value) {
  iterator position = find(tag);
  iterator limit = position;
  if (position != end())  ++limit;
  replace(position, limit, tag, value);
}

/// Print a header to the stream
/** Writes a header to an output stream as text in SAM format, @b without
a trailing newline character.
@relatesalso header */
std::ostream& operator<< (std::ostream& stream, const header& header);

/// Write the SAM representation of @a header to @a dest
/** @param dest  Character array to be written to; must have space for at least
header::sam_length() characters.
@param header  The header to be formatted.
@return  A pointer to the first unused character position in @a dest.
@relatesalso header */
char* format_sam(char* dest, const header& header);

/// Print a header field to the stream
/** @relatesalso header::tagfield */
std::ostream& operator<< (std::ostream& stream, const header::tagfield& field);

/// Print an iterator or const_iterator to the stream
std::ostream& operator<< (std::ostream& stream, header::const_iterator it);


/** @class sam::refsequence sam/header.h
    @brief Reference sequence record, corresponding to a single '@@SQ' header */
class refsequence : public header {
public:
  refsequence(const std::string& name, coord_t length, int index);

  ~refsequence() { }

  int index() const { return index_; }

  /** @name Reference sequence fields
  Accessors and modifiers for the defined reference sequence fields.  */
  //@{
  std::string name() const { return name_; }
  const char* name_c_str() const { return name_.c_str(); }
  coord_t length() const { return field<coord_t>("LN"); }
  std::string species() const { return field<std::string>("SP"); }
  std::string assembly() const { return field<std::string>("AS"); }
  std::string uri() const { return field<std::string>("UR"); }
  std::string checksum() const { return field<std::string>("M5"); }

  void set_name(const std::string& name) { set_field("SN", name); }
  void set_length(coord_t length) { set_field("LN", length); }
  void set_species(const std::string& species) { set_field("SP", species); }
  void set_assembly(const std::string& assembly) { set_field("AS", assembly); }
  void set_uri(const std::string& uri) { set_field("UR", uri); }
  void set_checksum(const std::string& sum) { set_field("M5", sum); }
  //@}

protected:
  virtual void sync();

private:
  friend class bamio;
  static std::string name_length_string(const std::string&, coord_t);

  friend class collection;
  refsequence(const std::string& nul_delimited_text, int index);

  std::string name_;
  int index_;
};


/** @class sam::readgroup sam/header.h
    @brief Read group record, corresponding to a single '@@RG' header */
class readgroup : public header {
public:
  ~readgroup() { }

  /** @name Read group fields
  Accessors and modifiers for the defined read group fields.  */
  //@{
  std::string id() const { return id_; }
  const char* id_c_str() const { return id_.c_str(); }
  std::string sample() const { return field<std::string>("SM"); }
  std::string library() const { return field<std::string>("LB"); }
  std::string description() const { return field<std::string>("DS"); }
  std::string unit() const { return field<std::string>("PU"); }
  coord_t median_isize() const { return field<coord_t>("PI"); }
  //@}

protected:
  virtual void sync();

private:
  friend class collection;
  readgroup(const std::string& nul_delimited_text);

  std::string id_;
};


/** @class sam::collection sam/header.h
    @brief Header information for a collection of SAM/BAM records */
class collection {
public:
  /// Construct an empty collection
  collection();

  /// Construct a copy of a collection, by copying all the headers within
  collection(const collection& collection);

  //  Destroy this collection object (not interesting enough to warrant ///)
  ~collection();

  /// Copy a collection, by copying all the headers within
  collection& operator= (const collection& collection);

  /** @name Container functionality
  Collections provide limited container-style access to their headers and
  reference sequences.

  The various @c iterator classes are <b>random access iterators</b> providing
  all the usual iterator functionality.  Dereferencing an @c iterator or
  @c const_iterator produces a sam::header, while dereferencing a
  @c ref_iterator or @c const_ref_iterator produces a sam::refsequence.  */
  //@{

  // @cond infrastructure
  template <bool condition, typename TrueType, typename FalseType>
  struct if_c { typedef TrueType type; };

  template <typename TrueType, typename FalseType>
  struct if_c<false, TrueType, FalseType> { typedef FalseType type; };

  template <typename ValueType, typename ContainerType, bool is_const = false>
  class ptr_container_iterator {
  public:
    // FIXME should be random_access_iterator_tag
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef ValueType value_type;
    typedef ptrdiff_t difference_type;
    typedef typename if_c<is_const, const ValueType*, ValueType*>::type pointer;
    typedef typename if_c<is_const, const ValueType&, ValueType&>::type
      reference;

    // This is either a copy constructor (for a mutable iterator, when is_const
    // is false) or a constructor (for a const_iterator, when is_const is true)
    // converting from the corresponding iterator.  In the latter case, the
    // compiler will generate a const_iterator copy constructor itself.
    ptr_container_iterator(const ptr_container_iterator<ValueType,
				   ContainerType, false>& it) : pit(it.pit) { }

    ptr_container_iterator() { }
    ~ptr_container_iterator() { }
    ptr_container_iterator& operator= (ptr_container_iterator it)
      { pit = it.pit; return *this; }

    // FIXME something about friend functions for these...
    bool operator== (ptr_container_iterator it) const { return pit == it.pit; }
    bool operator!= (ptr_container_iterator it) const { return pit != it.pit; }

    pointer operator-> () const { return *pit; }
    reference operator* () const { return **pit; }

    ptr_container_iterator& operator++ () { ++pit; return *this; }
    ptr_container_iterator operator++ (int)
      { ptr_container_iterator orig = *this; ++pit; return orig; }

    ptr_container_iterator& operator-- () { --pit; return *this; }
    ptr_container_iterator operator-- (int)
      { ptr_container_iterator orig = *this; --pit; return orig; }

  private:
    typedef typename
      if_c<is_const, typename ContainerType::const_iterator,
		     typename ContainerType::iterator>::type underlying_type;

    friend class collection;
    explicit ptr_container_iterator(const underlying_type& p) : pit(p) { }

    underlying_type pit;
  };

  typedef ptr_container_iterator<header, std::vector<header*> > iterator;
  typedef ptr_container_iterator<header, std::vector<header*>, true>
    const_iterator;
  typedef iterator::reference reference;
  typedef const_iterator::reference const_reference;
  typedef const_iterator::difference_type difference_type;

  typedef ptr_container_iterator<refsequence, std::vector<refsequence*> >
    ref_iterator;
  typedef ptr_container_iterator<refsequence, std::vector<refsequence*>, true>
    const_ref_iterator;
  // @endcond

  iterator begin() { return iterator(headers.begin()); }
  const_iterator begin() const { return const_iterator(headers.begin()); }

  iterator end() { return iterator(headers.end()); }
  const_iterator end() const { return const_iterator(headers.end()); }

  void push_back(const std::string& header_line);

  size_t size() const { return headers.size(); }
  bool empty() const { return headers.empty(); }
  void clear();

  ref_iterator ref_begin() { return ref_iterator(refseqs.begin()); }
  const_ref_iterator ref_begin() const
    { return const_ref_iterator(refseqs.begin()); }

  ref_iterator ref_end() { return ref_iterator(refseqs.end()); }
  const_ref_iterator ref_end() const
    { return const_ref_iterator(refseqs.end()); }

  size_t ref_size() const { return refseqs.size(); }
  bool ref_empty() const { return refseqs.empty(); }
  //@}

  /// @name Searching
  //@{
  refsequence& findseq(const std::string& name) { return findseq_(name); }
  const refsequence& findseq(const std::string& name) const
    { return findseq_(name); }
  refsequence& findseq(const char* name) { return findseq_(name); }
  /// Find reference sequence by @a name, which may be "*"
  const refsequence& findseq(const char* name) const { return findseq_(name); }

  refsequence& findseq(int index) { return findseq_(index); }
  /// Find reference sequence by @a index, where -1 <= @a index < ref_size()
  const refsequence& findseq(int index) const { return findseq_(index); }

  readgroup& findgroup(const std::string& id) { return findgroup_(id); }
  const readgroup& findgroup(const std::string& id) const
    { return findgroup_(id); }
  readgroup& findgroup(const char* id) { return findgroup_(id); }
  /// Find read group by @a id
  const readgroup& findgroup(const char* id) const { return findgroup_(id); }
  //@}

  // FIXME prob not public
  static collection& find(unsigned cindex) { return *collections[cindex]; }

private:
  friend class sambamio;
  friend class bamio;
  friend class samio;
  // @cond private
  friend std::ostream& operator<< (std::ostream&, const collection&);
  // @endcond

  void allocate_cindex();
  void free_cindex() { collections[cindex] = NULL; cindex = 0; }
  void reallocate_cindex() { collections[cindex] = NULL; allocate_cindex(); }

  void push_back(const std::string& nul_delimited_text, int flags);

  refsequence& findseq_(const std::string& name) const;
  refsequence& findseq_(const char* name) const;
  refsequence& findseq_(int index) const;
  readgroup& findgroup_(const std::string& id) const;

  typedef std::map<std::string, refsequence*> refname_map;
  typedef std::map<std::string, readgroup*> readgroup_map;

  int cindex;
  std::vector<header*> headers;
  std::vector<refsequence*> refseqs;
  refname_map refnames;
  bool refseqs_in_headers;
  readgroup_map rgroups;

  static std::vector<collection*> collections;
};

/// Print a collection of headers to the stream
/** Writes the headers to an output stream as text in SAM format, with each
header terminated by a newline character.  If @a stream has @c showpoint set,
also writes a representation of the reference list and other internal data.
@relatesalso collection */
std::ostream& operator<< (std::ostream& stream, const collection& headers);

} // namespace sam

#endif
