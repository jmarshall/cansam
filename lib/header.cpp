#include "sam/header.h"

#include <string>

#include "sam/exception.h"
#include "lib/utilities.h"

using std::string;

namespace sam {

// Infrastructure
// ==============

string header::type() const {
  if (! (str_.length() >= 3 && str_[0] == '@'))
    throw bad_format("Malformatted header type");

  return str_.substr(1, 2);
}

string header::tagfield::tag() const {
  if (! (tag_[0] && tag_[1] && colon_ == ':'))
    throw bad_format("Malformatted header field");

  return string(tag_, sizeof tag_);
}

string header::tagfield::value_str() const {
  const char* limit = nexttab(tag_);
  if (! (limit >= data_ && colon_ == ':'))
    throw bad_format("Malformatted header field");

  return string(data_, limit - data_);
}

template<> int header::tagfield::value<int>() const {
}

template<> coord_t header::tagfield::value<coord_t>() const {
}

size_t header::find_or_eos(const char* tag) const {
  char key[] = { '\t', tag[0], tag[1], ':' };
  size_t pos = str_.find(key, 3, sizeof key);
  return (pos != string::npos)? pos : str_.length();
}

header::const_iterator header::find_or_throw(const char* tag) const {
  size_t pos = find_or_eos(tag);
  if (pos == str_.length())
    throw bad_format(make_string()
	<< "Header field '" << tag[0] << tag[1] << "' not found");

  return const_iterator(cstr_ + pos);
}

header::iterator header::erase(iterator start, iterator limit) {
  size_t pos = start.ptr - cstr_;
  str_.erase(pos, limit.ptr - start.ptr);
  sync();
  return iterator(cstr_ + pos);
}

void header::clear() {
  str_.erase(3);
  sync();
}

int header::erase(const char* key_tag) {
  int count = 0;
  iterator it = begin();
  while (it != end())
    if (it->tag_equals(key_tag))  it = erase(it), count++;
    else  ++it;

  return count;
}

header::iterator
header::replace_string(size_t pos, size_t length,
		       const char* tag, const char* value, int value_length) {
  char key[] = { '\t', tag[0], tag[1], ':' };

  if (length >= sizeof key) {
    str_.replace(pos, sizeof key, key, sizeof key);
    str_.replace(pos + sizeof key, length - sizeof key, value, value_length);
  }
  else if (pos == str_.length()) {
    str_.reserve(str_.length() + sizeof key + value_length);
    str_.append(key, sizeof key);
    str_.append(value, value_length);
  }
  else {
    str_.insert(pos, sizeof key + value_length - length, '#');
    str_.replace(pos, sizeof key, key, sizeof key);
    str_.replace(pos + sizeof key, value_length, value, value_length);
  }

  sync();
  return iterator(cstr_ + pos);
}

header::iterator
header::replace_(size_t pos, size_t length, const char* tag, int value) {
  char buffer[format::int_digits];
  char* buflim = format::decimal(buffer, value);
  return replace_string(pos, length, tag, buffer, buflim - buffer);
}

#if 0
// FIXME Need to think about whether replace(const_iterator) should be able
// to update the destination's tag.  Ditto for alignments.
header::iterator
header::replace_(size_t pos, size_t length,
		 const char* tag, const_iterator value) {
  const_iterator next = value;
  ++next;
  str_.replace(pos, length, value.ptr, next.ptr - value.ptr);

  sync();
  return iterator(cstr_ + pos);
}
#endif


// Reference sequences
// ===================

string refsequence::name_length_string(const string& name, coord_t length) {
  string s;
  // FIXME Should be format::coord_digits (and similarly below)
  s.reserve(7 + name.length() + 4 + format::int_digits);
  s += "@SQ\tSN:";
  s += name;
  s += "\tLN:";
  char buffer[format::int_digits];
  s.append(buffer, format::decimal(buffer, length));

  return s;
}

refsequence::refsequence(const string& name, coord_t length, int index)
  : header(name_length_string(name, length)),
    name_(name), index_(index), visible_(false) {
}

void refsequence::sync() {
  header::sync();
  name_ = field<string>("SN");

  // FIXME if (this == &unmapped_refseq) throw barf;

  // If a reference has ever been modified, it should appear in the header.
  visible_ = true;
}

} // namespace sam
