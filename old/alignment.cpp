#define __STDC_LIMIT_MACROS
#include "sam/alignment.h"

#include <string>

#if !(defined UINT8_MAX && defined UINT16_MAX && defined INT32_MAX)
#error This implementation requires 8, 16, and 32-bit integer types
#endif

namespace sam {

// For order(); 0, FIRST_IN_PAIR, SECOND_IN_PAIR, FIRST_IN_PAIR|SECOND_IN_PAIR.
const int alignment::order_value[4] = { 0, -1, +1, 0 };

#if 0
std::string alignment::qname() const {
  // FIXME if (qname_.empty()) { stuff; }
  return qname_;
}

std::string& alignment::qname() {
#if 0
  if (qname_.empty())
    qname_.assign(&buffer[36/*READNAME_offset*/], cigar_offset - READNAME_offset - 1);
#endif
  return qname_;
}

//============
coord_t alignment::pos() const {
  return pos_;
}

std::string alignment::rname() const {
  return rname_;
}

std::string alignment::cigar() const {
  return cigar_;
}

std::string alignment::mate_rname() const {
  return mate_rname_;
}

coord_t alignment::mate_pos() const {
  return mate_pos_;
}

std::string alignment::seq() const {
  return seq_;
}

std::string alignment::qual() const {
  return qual_;
}

const tagfield::array& alignment::auxen() const {
  return auxen_;
}

//============

// FIXME or should it take (begin, length)?
int alignment::to_flags(const char* s, const char* lim) {
  size_t len = lim - s;

  int val = 0;

  if (len == 0) {
    // The field is empty, which is mostly invalid and unlikely,
    // so returning zero is not much worse than anything else.
  }
  else if (s[0] == '0' && len >= 2 && (s[1] == 'x' || s[1] == 'X')) {
    for (size_t i = 2; i < len; i++, s++) {
      val *= 16;
      if (*s >= 'a')  val += *s - 'a' + 10;
      else if (*s >= 'A')  val += *s - 'A' + 10;
      else  val += *s - '0';
    }
  }
  else if (s[0] >= '0' && s[1] <= '9') {
    int base = (s[0] == '0')? 8 : 10;
    for (size_t i = 0; i < len; i++, s++) {
      val *= base;
      val += *s - '0';
    }
  }
  else {
    // FIXME One day there'll be a defined text representation
  }

  return val;
}

int alignment::to_flags(const std::string& str) {
  const char* s = str.data();
  return to_flags(s, s + str.length());
}


// FIXME move us
inline bool tag_eq(const char* s1, const char* s2) {
  return s1[0] == s2[0] && s1[1] == s2[1];
}

tagfield::array::const_iterator tagfield::array::find(const char* tag) const {
  for (const_iterator it = begin(); it != end(); ++it)
    if (tag_eq(it->tag, tag))
      return it;

  return end();
}

tagfield::array::iterator tagfield::array::find(const char* tag) {
  for (iterator it = begin(); it != end(); ++it)
    if (tag_eq(it->tag, tag))
      return it;

  return end();
}
// end move us


char alignment::has(const char* tag) const {
  tagfield::array::const_iterator it = auxen_.find(tag);
  return (it != auxen_.end())? it->type : '\0';
}

std::string alignment::aux(const char* tag, const std::string& defval) const {
  tagfield::array::const_iterator it = auxen_.find(tag);
  return (it != auxen_.end())? it->str : defval;
}

std::string& alignment::aux(const char* tag) {
  tagfield::array::iterator it = auxen_.find(tag);
  if (it != auxen_.end())
    return it->str;

  auxen_.push_back(tagfield(tag, 'Z'));
  return auxen_.back().str;
}

int alignment::aux_int(const char* tag, int default_value) const {
  tagfield::array::const_iterator it = auxen_.find(tag);
  return (it != auxen_.end())? it->i : default_value;
}

int& alignment::aux_int(const char* tag) {
  tagfield::array::iterator it = auxen_.find(tag);
  if (it != auxen_.end())
    return it->i;

  auxen_.push_back(tagfield(tag, 'i'));
  return auxen_.back().i;
}
#endif

} // namespace sam
