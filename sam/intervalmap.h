/// @file sam/intervalmap.h
/// Classes for sequence intervals and interval containers

/*  Copyright (C) 2011 Genome Research Ltd.

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

#ifndef CANSAM_INTERVALMAP_H
#define CANSAM_INTERVALMAP_H

#include <functional>
#include <iosfwd>
#include <iterator>
#include <map>
#include <string>
#include <utility>

#include <stdint.h>

#include <iostream> // FIXME NUKEME

#include "sam/types.h"

namespace sam {

/** @class sam::interval sam/intervalmap.h
    @brief Interval within an unspecified sequence
*/
class interval {
public:
  /// Construct an empty interval
  interval() : zstart_(0), zlimit_(0) { }

  /// Construct a zero-based, half-open interval
  interval(coord_t zstart, coord_t end) : zstart_(zstart), zlimit_(end) { }

  /// Construct a copy of an interval
  interval(const interval& i) : zstart_(i.zstart_), zlimit_(i.zlimit_) { }

  //  Destroy this interval object (not interesting enough to warrant ///)
  ~interval() { }

  /// Copy an interval
  interval& operator= (const interval& i)
    { zstart_ = i.zstart_; zlimit_ = i.zlimit_; return *this; }

  /** @name Field accessors
  Two variants are provided for each field: @c %start() et al return 1-based
  coordinates, while @c %zstart() et al return the same positions in 0-based
  coordinates.  */
  //@{
  coord_t start() const  { return zstart_+1; }
  coord_t zstart() const { return zstart_; }

  coord_t end() const  { return zlimit_; }
  coord_t zend() const { return zlimit_-1; }

  coord_t limit() const  { return zlimit_+1; }
  coord_t zlimit() const { return zlimit_; }

  coord_t length() const { return zlimit_ - zstart_;  }
  //@}

  /// @name Field modifiers
  //@{
  void set_start(coord_t start)   { zstart_ = start-1; }
  void set_zstart(coord_t zstart) { zstart_ = zstart; }

  void set_end(coord_t end)   { zlimit_ = end; }
  void set_zend(coord_t zend) { zlimit_ = zend+1; }

  void set_limit(coord_t limit)   { zlimit_ = limit-1; }
  void set_zlimit(coord_t zlimit) { zlimit_ = zlimit; }
  //@}

protected:
  // @cond infrastructure
  // Represented as a zero-based half-open interval [zstart_, zlimit_).
  int32_t zstart_;
  int32_t zlimit_;
  // @endcond
};

/// Compare intervals
/** The intervals are ordered by starting coordinate.
@relatesalso interval */
inline bool operator< (const interval& a, const interval& b)
  { return a.zstart() < b.zstart(); }

/// Compare intervals, similarly to @c operator< above
/** @relatesalso interval */
inline bool operator> (const interval& a, const interval& b) { return b < a; }

// TODO  Perhaps pick an operator for this
inline bool overlaps(const interval& a, const interval& b)
  { return a.zstart() < b.zlimit() && b.zstart() < a.zlimit(); }

/** @class sam::seqinterval sam/intervalmap.h
    @brief Interval within a named sequence
*/
class seqinterval : public interval {
public:
  /// Construct an empty seqinterval
  seqinterval() { }

  /// Construct a zero-based, half-open seqinterval
  seqinterval(const std::string& name, coord_t zstart, coord_t end)
    : interval(zstart, end), name_(name) { }

  /// Construct a zero-based, half-open seqinterval
  seqinterval(const char* name, coord_t zstart, coord_t end)
    : interval(zstart, end), name_(name) { }

  /// Construct a copy of a seqinterval
  seqinterval(const seqinterval& i) : interval(i), name_(i.name_) { }

  //  Destroy this seqinterval object (not interesting enough to warrant ///)
  ~seqinterval() { }

  /// Copy a seqinterval
  seqinterval& operator= (const seqinterval& i)
    { name_ = i.name_; zstart_ = i.zstart_; zlimit_ = i.zlimit_; return *this; }

  /// @name Field accessors
  //@{
  std::string name() const { return name_; }
  const char* name_c_str() const { return name_.c_str(); }
  //@}

  /// @name Field modifiers
  //@{
  void set_name(const std::string& name) { name_ = name; }
  //@}

private:
  std::string name_;
};

/// Print an interval to the stream
/** Writes an interval to a stream in the format "START-END", with 1-based
coordinates.
@relatesalso interval */
std::ostream& operator<< (std::ostream& stream, const interval& i);

/// Print a seqinterval to the stream
/** Writes a seqinterval to a stream in the format "NAME:START-END", with
1-based coordinates.
@relatesalso seqinterval */
std::ostream& operator<< (std::ostream& stream, const seqinterval& i);


// @cond infrastructure
class interval_tree_base {
public:
  class node {
  private:
    node* parent;
    node* left;
    node* right;
    int32_t max_zlimit;
    enum { red, black} colour;
  public:
    const interval first;

  protected:
    node(const interval& i) : parent(&nil), left(&nil), right(&nil),
			      max_zlimit(i.zlimit()), colour(red), first(i) { }
    ~node() { }

  private:
    friend class interval_tree_base;
    template<typename T> friend class interval_tree;  // for delete_tree() only
    node();
  };

  static bool is_red(const node* x) { return x->colour == node::red; }
  static bool is_left_child(const node* x)  { return x == x->parent->left; }
  static bool is_right_child(const node* x) { return x == x->parent->right; }

  void rotate_left(node* x);
  void rotate_right(node* x);

  static node* minimum(node* x)
    { while (x->left != &nil)  x = x->left;  return x; }

  static node* minimum_perhaps_intersecting(node* x, const interval& i) {
    while (i.zstart() < x->left->max_zlimit)  x = x->left;
    return x;
  }

  static node* successor_perhaps_intersecting(node* x, const interval& i) {
    if (x->right != &nil && x->first.zstart() < i.zlimit())
      return minimum_perhaps_intersecting(x->right, i);
    else {
      while (x->parent != &nil && is_right_child(x))  x = x->parent;
      return x->parent;
    }
  }

  static node* first_intersecting(node* x, const interval& i);
  static node* next_intersecting(node* x, const interval& i);

  void dump(const char* message) const {
    std::clog << "Interval tree " << this << " (" << message << ")\n";
    dump(0, 'T', root, &nil);
    std::clog << "---- end interval tree\n";
  }

  void dump_intersecting_r(const interval& i)
    { if (root != &nil)  dump_intersecting_r(root, i); }

  void dump_intersecting_i(const interval& i)
    { if (root != &nil)  dump_intersecting_i(root, i); }

protected:
  interval_tree_base() : root(&nil) { }
  ~interval_tree_base() { }
  node* insert(node* z);

  node* root;

private:
  static node nil;

  void dump(int level, char side, const node* p, const node* pparent) const;
  void dump_intersecting_r(node* x, const interval& i);
  void dump_intersecting_i(node* x, const interval& i);
};

template <typename MappedType>
class interval_tree : public interval_tree_base {
public:
  interval_tree() { }
  ~interval_tree() { delete_tree(root); }

private:
  struct pair_node : public node {
    pair_node(const interval& i, const MappedType& v) : node(i), second(v) { }

    MappedType second;
  };

  void delete_tree(node* x);

public:
  typedef pair_node value_type;

  class iterator : public std::iterator<std::forward_iterator_tag, value_type> {
  public:
    iterator(node* x) : ptr(x) { }
    iterator(node* x, const interval& i) : ptr(x), key(i) { }
    iterator(const iterator& it) : ptr(it.ptr), key(it.key) { }
    ~iterator() { }
    iterator& operator= (iterator it)
      { ptr = it.ptr; key = it.key; return *this; }

    bool operator== (iterator it) const { return ptr == it.ptr; }
    bool operator!= (iterator it) const { return ptr != it.ptr; }

    value_type& operator* () const  { return *static_cast<pair_node*>(ptr); }
    value_type* operator-> () const { return  static_cast<pair_node*>(ptr); }

    iterator& operator++ () { ptr = next_intersecting(ptr, key); return *this; }
    iterator operator++ (int)
      { node* orig = ptr;  ptr = next_intersecting(ptr, key);
	return iterator(orig, key); }

  private:
    node* ptr;
    interval key;
  };

  iterator begin() { return iterator(minimum()); }
  iterator end() { return iterator(&nil); }

  iterator insert(const interval& i, const MappedType& v)
    { return iterator(interval_tree_base::insert(new pair_node(i, v))); }

  std::pair<iterator, iterator> intersecting_range(const interval& i)
    { return make_pair(iterator(first_intersecting(root, i), i), end()); }
};

template <typename MappedType>
void interval_tree<MappedType>::delete_tree(node* x) {
  while (x != &nil) {
    if (x->right != &nil)  delete_tree(x->right);
    node* y = x->left;
    delete static_cast<pair_node*>(x);
    x = y;
  }
}
// @endcond

/** @class sam::interval_multimap sam/intervalmap.h
    @brief Associative container keyed by sequence intervals
*/
template <typename MappedType>
class interval_multimap {
public:
  // @cond infrastructure
  typedef MappedType mapped_type;
  typedef typename interval_tree<MappedType>::value_type value_type;
  typedef typename interval_tree<MappedType>::iterator iterator;
  typedef std::pair<iterator, iterator> iterator_pair;
  // @endcond

  /// Construct an empty container
  interval_multimap() { }

  //  Destroy this container object (not interesting enough to warrant ///)
  ~interval_multimap() { }

  /// Insert (a copy of) @a value
  iterator insert(const std::pair<seqinterval, MappedType>& value)
    { return trees[value.first.name()].insert(value.first, value.second); }

  /// Returns a range [@em first, @em last) of intervals intersecting @a i
  iterator_pair intersecting_range(const seqinterval& i)
    { return trees[i.name()].intersecting_range(i); }

  // FIXME NUKE-ME
  void dump_intersecting(const seqinterval& i) {
    std::cout << i << ':'; trees[i.name()].dump_intersecting_r(i); std::cout << '\n';
    std::cout << i << ':'; trees[i.name()].dump_intersecting_i(i); std::cout << '\n';
    }

private:
  // Prevent copy construction and assignment
  interval_multimap(const interval_multimap&);
  interval_multimap& operator= (const interval_multimap&);

  std::map<std::string, interval_tree<MappedType> > trees;
};

} // namespace sam

#endif
