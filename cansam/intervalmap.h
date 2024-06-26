/// @file cansam/intervalmap.h
/// Classes for sequence interval containers

/*  Copyright (C) 2011-2013 Genome Research Ltd.

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

#include <iterator>
#include <map>
#include <string>
#include <utility>

#include <stdint.h>

#include <iostream> // FIXME NUKEME

#include "cansam/interval.h"

namespace sam {

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

  static node nil;

private:
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

  class iterator {
  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef pair_node value_type;
    typedef std::ptrdiff_t difference_type;
    typedef pair_node* pointer;
    typedef pair_node& reference;

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

  iterator begin() { return iterator(minimum(root)); }
  iterator end() { return iterator(&nil); }

  iterator insert(const interval& i, const MappedType& v)
    { return iterator(interval_tree_base::insert(new pair_node(i, v))); }

  std::pair<iterator, iterator> intersecting_range(const interval& i)
    { return std::make_pair(iterator(first_intersecting(root, i), i), end()); }
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

/** @class sam::interval_multimap cansam/intervalmap.h
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

  /// Returns whether the container is empty
  bool empty() const { return trees.empty(); }

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
