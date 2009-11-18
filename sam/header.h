/** @file   sam/header.h
    @brief  Classes for SAM/BAM headers
*/

#ifndef CANSAM_HEADER_H
#define CANSAM_HEADER_H

#include <iosfwd>
#include <string>
#include <vector>

namespace sam {

/** @class sam::header_field sam/header.h
    @brief Header field with its two-character tag
*/
struct header_field {
  header_field(const std::string& tag0, const std::string& value0)
    : tag(tag0), value(value0) { }

  std::string tag;
  std::string value;
};

/** @class sam::header sam/header.h
    @brief SAM/BAM %header record, representing a single '@@' %header line

This class represents the %header fields as a vector of header_field items.
This is by inheritance, so a header object itself is amenable to all the
usual vector operations, such as iterating through the fields with @c begin()
and @c end(), and modifying them with @c erase(), @c %push_back(), etc.

A SAM %header probably ought not to have two fields with the same tag,
though this is not explicitly stated in the specification and is not enforced
by this class.  Thus find() and field() in fact return the @b first field with
the specified tag.  In the unlikely event that you expect a tag to occur in
several fields, you should iterate until @c end() checking for other instances
of the tag yourself.  */
class header : public std::vector<header_field> {
public:
  /// Construct a %header with an empty type code and no fields.
  header() { }

  /// Construct a %header by splitting up a tab-separated text string.
  explicit header(const std::string& line) { assign(line); }

  /// Destroy the %header object.
  ~header() { }

  /// Assign to a %header by splitting up a tab-separated text string.
  header& assign(const std::string& line);

  /// Returns the tab-separated string representing the %header.
  std::string str() const;

  /// Returns the %header's two-character type code.
  std::string type() const { return type_; }

  /// Returns the field with the specified @a tag, or @c end() if not found.
  const_iterator find(const char* tag) const;

  /// Returns the field with the specified @a tag, or @c end() if not found.
  iterator find(const char* tag);

  /** @brief Returns the text of the field with the specified @a tag,
	     or @a default_value if not found.  */
  std::string field(const char* tag,
		    const std::string& default_value = std::string()) const;

  // FIXME Maybe also field_int() etc
  
  using std::vector<header_field>::push_back;
  /** @brief [In addition to the inherited @c %push_back]
	     Add a new field at the end of the %header.  */
  void push_back(const char* tag, const std::string& value)
    { push_back(header_field(tag, value)); }

private:
  std::string type_;
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

} // namespace sam

#endif
