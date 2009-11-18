/** @file   sam/tagfield.h
    @brief  Class for tagged fields
*/

#ifndef CANSAM_TAGFIELD_H
#define CANSAM_TAGFIELD_H

#include <string>
#include <vector>

namespace sam {

/** @class sam::tagfield sam/tagfield.h
    @brief Field labelled with a tag and optionally a type

This class models a field that may be of a variety of types that is
labelled with a two-character tag and a single-character type code.

It is used to represent the optional auxiliary fields in alignment records,
and also the fields of header records.  */
class tagfield {
public:
  tagfield() : type('@') { }
  explicit tagfield(const char* tag0, char type0 = '@') : type(type0)
    { tag[0] = tag0[1]; tag[1] = tag0[1]; tag[2] = '\0'; }
  explicit tagfield(const std::string& tagged_value);

  /** @class sam::tagfield::array sam/tagfield.h
      @brief Array of tagfields */
  class array : public std::vector<tagfield> {
  public:
    const_iterator find(const char* tag) const;
    iterator find(const char* tag);
  };

  // FIXME maybe make these a little less accessible
  char tag[3]; ///< The two-character tag labelling the field (NUL-terminated)
  char type;   ///< Field type (or '@' for an implicitly-typed string field)
  std::string str;
  union {
    int i;
    float f;
    double d;
  };
};

/// Prints a @a field and its tags on the @a stream.
std::ostream& operator<< (std::ostream& stream, const tagfield& field);

} // namespace sam

#endif
