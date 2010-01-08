/** @file   sam/header.h
    @brief  Classes for SAM/BAM headers
*/

#ifndef CANSAM_HEADER_H
#define CANSAM_HEADER_H

#include <iosfwd>
#include <string>
#include <vector>

#include "sam/types.h"

namespace sam {

/** @class sam::header sam/header.h
    @brief SAM/BAM header record, representing a single '@@' header line

FIXME Doco update, it's not a std::vector any more!

This class represents the header fields as a vector of header_field items.
This is by inheritance, so a header object itself is amenable to all the
usual vector operations, such as iterating through the fields with @c begin()
and @c end(), and modifying them with @c erase(), @c %push_back(), etc.

A SAM header probably ought not to have two fields with the same tag,
though this is not explicitly stated in the specification and is not enforced
by this class.  Thus find() and field() in fact return the @b first field with
the specified tag.  In the unlikely event that you expect a tag to occur in
several fields, you should iterate until @c end() checking for other instances
of the tag yourself.  */
class header {
private:
  // @cond private
  class header_field {
  public:
    header_field(const char* tag, const std::string& value0)
      : value(value0) { tag_[0] = tag[0]; tag_[1] = tag[1]; tag_[2] = '\0'; }

    std::string tag() const { return tag_; }
    std::string field() const { return value; }
    int field_int() const { return 37; } // FIXME

    void set_field(const std::string& newvalue);
    void set_field(int newvalue);

  private:
    char tag_[3];
    std::string value;
    };

  typedef std::vector<header_field> field_array;
  // @endcond

public:
  // @cond private
  typedef field_array::size_type size_type;
  typedef field_array::const_iterator const_iterator;
  typedef field_array::iterator iterator;
  // @endcond

  /// Construct a header with an empty type code and no fields
  header() { }

  /// Construct a header by splitting up a tab-separated text string
  explicit header(const std::string& line) { assign(line); }

  /// @brief Allocate a new particular kind of header by splitting up
  /// a tab-separated text string
  static header* new_header(const std::string& line);

  /// Destroy this header object
  virtual ~header() { }

  /// Assign to a header by splitting up a tab-separated text string
  header& assign(const std::string& line);

  /// Returns the tab-separated string representing the header
  std::string str() const;

  /// Returns the header's two-character type code
  std::string type() const { return type_; }

  /// @brief Text of the field with the given @a tag, or @a default_value
  /// if there is no such field
  std::string field(const char* tag,
		    const std::string& default_value = std::string()) const;

  /// @brief The field with the given @a tag as an integer, or @a default_value
  /// if there is no such field
  int field_int(const char* tag, int default_value = 0) const;


  void set_type(const std::string& type) { type_ = type; sync(); }

  void set_field(const char* tag, const std::string& value);
  void set_field(const char* tag, int value);

  void set_field(iterator pos, const std::string& value);
  void set_field(iterator pos, int value);

  /** @name Collection Public Members
  Headers provide limited collection-style access to their fields,
  similar to @c std::vector.  */
  //@{
  /// Returns whether the header has no fields at all
  bool empty() const { return fields.empty(); }

  /// Returns the number of fields contained in the header
  size_type size() const { return fields.size(); }

  const_iterator begin() const { return fields.begin(); }
  iterator begin() { return fields.begin(); }

  const_iterator end() const { return fields.begin(); }
  iterator end() { return fields.begin(); }

  /// Returns the field with the specified @a tag, or @c end() if not found.
  const_iterator find(const char* tag) const;

  /// Returns the field with the specified @a tag, or @c end() if not found.
  iterator find(const char* tag);

  /// Erase all fields from the header
  void clear() { fields.clear(); sync(); }

  iterator erase(iterator pos) { pos = fields.erase(pos); sync(); return pos; }
  void erase(const char* tag);

  iterator insert(iterator it, const char* tag, const std::string& value);
  iterator insert(iterator it, const char* tag, int value);

  void push_back(const char* tag, const std::string& value)
    { fields.push_back(header_field(tag, value)); sync(); }

  void push_back(const char* tag, int value);
  //@}

private:
  virtual void sync() { }

  std::string type_;
  std::vector<header_field> fields;
};

/** @brief Returns the kind (header/alignment) of the next available record
@param stream  The input stream to be peeked at.
@return  Either '@' for a header, 'A' for an alignment record, or EOF.
@relatesalso header  */
int peek(std::istream& stream);

/** @brief Reads a %header from the stream.
@details Extracts a single SAM-formatted @c @@%header line from the input
@a stream into @a hdr.  Sets @c failbit if the first line available does
not start with an '@' character, i.e., if it is not a %header line.
@relatesalso header  */
std::istream& operator>> (std::istream& stream, header& hdr);

/** @brief Prints a %header on the stream.
@details   Writes @a hdr to the output @a stream as text in SAM format,
@b without a trailing newline character.
@relatesalso header  */
std::ostream& operator<< (std::ostream& stream, const header& hdr);


/** @class sam::reference sam/header.h
    @brief Reference sequence record, corresponding to a single '@@SQ' header */
class reference : public header {
public:
  reference() { }

  reference(const std::string& name, coord_t length);
//    : name_(name), length_(length)
//    { name_it = find("SN"); }

  ~reference() { }

  /// Name of the sequence
  std::string name() const { return name_it->field(); }

  /// Sequence length
  coord_t length() const { return field_int("LN"); }

  /// Species
  std::string species() const { return field("SP"); }

  /// Genome assembly
  std::string assembly() const { return field("AS"); }

  /// Sequence location URI
  std::string uri() const { return field("UR"); }

  /// MD5 checksum of the sequence text
  std::string checksum() const { return field("M5"); }

  void set_name(const std::string& name) { set_field("SN", name); }
  void set_length(coord_t length) { set_field("LN", length); }

private:
  virtual void sync();

  iterator name_it;
};

} // namespace sam

#endif
