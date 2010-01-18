/** @file   sam/header.h
    @brief  Classes for SAM/BAM headers
*/

#ifndef CANSAM_HEADER_H
#define CANSAM_HEADER_H

#include <string>
#include <iterator>
#include <iosfwd>

#include "sam/types.h"

namespace sam {

/** @class sam::header sam/header.h
    @brief SAM/BAM header record, representing a single '@@' header line

This class blah FIXME blah collection-style blah, so a header object itself
is amenable to the usual collection operations, such as iterating through the
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
  explicit header(const std::string& line) : str_(line), cstr_(str_.c_str()) { }
  virtual ~header() { }

  // FIXME Infrastructure: copy ctor, assign, swap, reserve

  /// The header's two-character type code
  std::string type() const;

  /// Returns whether this header's type code is the same as the
  /// given @a key_type
  bool type_equals(const char* key_type) const
    { return cstr_[0]=='@' && cstr_[1]==key_type[0] && cstr_[2]==key_type[1]; }

  std::string text() const;

  /// Returns the (tab-separated) string representing this header
  std::string str() const { return str_; }

  std::string field_str(const char* tag) const;

  template <typename ValueType>
  ValueType field(const char* tag) const
    { return find_or_throw(tag)->template value<ValueType>(); }

  template <typename ValueType>
  ValueType field(const char* tag, ValueType default_value) const;

  /** @name Fields as a collection
  Headers provide limited collection-style access to their fields.

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
    std::string value_str() const;

    template <typename ValueType>
    ValueType value() const;

    /// Returns whether the field's tag is the same as the given @a key_tag
    bool tag_equals(const char* key_tag) const
      { return tag_[0] == key_tag[0] && tag_[1] == key_tag[1]; }

  private:
    friend class const_iterator;
    friend std::ostream& operator<< (std::ostream&, const tagfield&);

    inline static const char* nexttab(const char* s)
      { char c; do { c = *++s; } while (c != '\0' && c != '\t');  return s; }

    inline static const char* prevtab(const char* s)
      { while (*--s != '\t') { }  return s; }

    char tab_;
    char tag_[2];
    char colon_;
    char data_[1];
  };

  // @cond infrastructure
  class const_iterator :
    public std::iterator<std::bidirectional_iterator_tag, tagfield,
			 ptrdiff_t, const tagfield*, const tagfield&> {
  public:
    const_iterator() { }
    const_iterator(const const_iterator& it) : ptr(it.ptr) { }
    const_iterator& operator= (const_iterator it) {ptr = it.ptr; return *this;}
    ~const_iterator() { }

    bool operator== (const_iterator it) const { return ptr == it.ptr; }
    bool operator!= (const_iterator it) const { return ptr != it.ptr; }

    const tagfield* operator-> () const
      { return reinterpret_cast<const tagfield*>(ptr); }
    const tagfield& operator* () const { return *operator->(); }

    const_iterator& operator++ () {ptr = tagfield::nexttab(ptr); return *this;}
    const_iterator operator++ (int)
      { const char* orig = ptr; ++(*this); return const_iterator(orig); }

    const_iterator& operator-- () {ptr = tagfield::prevtab(ptr); return *this;}
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
    iterator& operator= (iterator it) { ptr = it.ptr; return *this; }
    ~iterator() { }

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

  iterator find(const char* tag) { return iterator(cstr_ + find_(tag)); }
  const_iterator find(const char* tag) const
    { return const_iterator(cstr_ + find_(tag)); }

  bool empty() const { return str_.length() <= 3; }

  template <typename ValueType>
  void push_back(const char* tag, ValueType value)
    { replace_(str_.length(), 0, tag, value); }

  template <typename ValueType>
  iterator insert(iterator position, const char* tag, ValueType value)
    { return replace_(position.ptr - cstr_, 0, tag, value); }

  iterator erase(iterator position) { return erase(position, next(position)); }
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
    { return replace(position, next(position), NULL, value); }

  /// Erase all fields with the given @a tag
  /** @return The number of fields erased (usually no more than 1).  */
  int erase(const char* tag);

protected:
  // @cond private
  virtual void sync() { cstr_ = str_.c_str(); }
  // @endcond

private:
  size_t find_(const char* tag) const;
  const_iterator find_or_throw(const char* tag) const;

  iterator replace_string(size_t pos, size_t length,
			  const char* tag, const char* value, int value_length);

  iterator replace_(size_t pos, size_t length,
		    const char* tag, const std::string& value)
    { return replace_string(pos, length, tag, value.data(), value.length()); }

  iterator replace_(size_t pos, size_t length,
		    const char* tag, const char* value)
    { return replace_string(pos, length, tag,
			    value, std::char_traits<char>::length(value)); }

  iterator replace_(size_t pos, size_t length, const char* tag, int value);
  iterator replace_(size_t pos, size_t length,
		    const char* tag, const_iterator value);

  std::string str_;
  const char* cstr_;
};

template <> std::string header::tagfield::value<std::string>() const { return value_str(); }
template <> int header::tagfield::value<int>() const { return 4; }
template <> coord_t header::tagfield::value<coord_t>() const { return 5; }

template <typename ValueType>
ValueType header::field(const char* tag, ValueType default_value) const {
  const_iterator it = find(tag);
  return (it != end())? it->value<ValueType>() : default_value;
}

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
@relatesalso header  */
std::ostream& operator<< (std::ostream& stream, const header& hdr);

/// Print a header field to the stream
/** @relatesalso header::tagfield */
std::ostream& operator<< (std::ostream& stream, const header::tagfield& field);

/// Print an iterator or const_iterator to the stream
std::ostream& operator<< (std::ostream& stream, header::const_iterator it);


/** @class sam::reference sam/header.h
    @brief Reference sequence record, corresponding to a single '@@SQ' header */
class reference : public header {
public:
  reference(const std::string& name, coord_t length);
//    : name_(name), length_(length)
//    { name_it = find("SN"); }

  ~reference() { }

  /** @name Reference sequence fields
  Accessors and modifiers for the defined reference sequence fields.  */
  //@{
  std::string name() const { return name_; }
  const char* name_c_str() const { return name_.c_str(); }
  coord_t length() const { return field<coord_t>("LN"); }
  std::string species() const { return field_str("SP"); }
  std::string assembly() const { return field_str("AS"); }
  std::string uri() const { return field_str("UR"); }
  std::string checksum() const { return field_str("M5"); }

  void set_name(const std::string& name) { set_field("SN", name); }
  void set_length(coord_t length) { set_field("LN", length); }
  void set_species(const std::string& species) { set_field("SP", species); }
  void set_assembly(const std::string& assembly) { set_field("AS", assembly); }
  void set_uri(const std::string& uri) { set_field("UR", uri); }
  void set_checksum(const std::string& sum) { set_field("M5", sum); }
  //@}

private:
  virtual void sync();

  std::string name_;
  bool visible_;
};

} // namespace sam

#endif
