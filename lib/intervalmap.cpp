/*  intervalmap.cpp -- Classes for sequence intervals and interval containers.

    Copyright (C) 2011 Genome Research Ltd.

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

#include "sam/intervalmap.h"

#include <ostream>

#include <iostream>  // FIXME NUKE-ME

namespace sam {

namespace {

inline int32_t max(int32_t i, int32_t j, int32_t k) {
  int32_t m = i;
  if (j > m)  m = j;
  if (k > m)  m = k;
  return m;
}

} // unnamed namespace


std::ostream& operator<< (std::ostream& stream, const interval& i) {
  return stream << i.start() << '-' << i.end();
}

std::ostream& operator<< (std::ostream& stream, const seqinterval& i) {
  return stream << i.name() << ':' << i.start() << '-' << i.end();
}


// This constructor is used only for the nil singleton, which is set up
// such that various tests conveniently fail: it has no children, isn't red,
// and no zstart precedes its max_zlimit.
interval_tree_base::node::node()
  : left(&nil), right(&nil), max_zlimit(-1), colour(black) {
}

interval_tree_base::node interval_tree_base::nil;

void interval_tree_base::rotate_left(node* x) {
  node* y = x->right;

  x->right = y->left;
  if (y->left != &nil)  y->left->parent = x;

  y->parent = x->parent;
  if (x->parent == &nil)  root = y;
  else if (is_left_child(x))  x->parent->left = y;
  else  x->parent->right = y;

  y->left = x;
  x->parent = y;

  y->max_zlimit = x->max_zlimit;
  x->max_zlimit = max(x->left->max_zlimit, x->first.zlimit(),
		      x->right->max_zlimit);
}

void interval_tree_base::rotate_right(node* y) {
  node* x = y->left;

  y->left = x->right;
  if (x->right != &nil)  x->right->parent = y;

  x->parent = y->parent;
  if (y->parent == &nil)  root = x;
  else if (is_left_child(y))  y->parent->left = x;
  else  y->parent->right = x;

  x->right = y;
  y->parent = x;

  x->max_zlimit = y->max_zlimit;
  y->max_zlimit = max(y->left->max_zlimit, y->first.zlimit(),
		      y->right->max_zlimit);
}

interval_tree_base::node* interval_tree_base::insert(node* z) {
  node* y = &nil;
  node* x = root;
  while (x != &nil) {
    y = x;
    x = (z->first < x->first)? x->left : x->right;
  }

  z->parent = y;
  if (y == &nil)  root = z;
  else if (z->first < y->first)  y->left = z;
  else  y->right = z;

  // The new node has been inserted; now fix up the max property...
  while (y != &nil && y->max_zlimit < z->max_zlimit) {
    y->max_zlimit = z->max_zlimit;
    y = y->parent;
  }

  // ...and fix up the red-black properties.
  while (is_red(z->parent))
    if (is_left_child(z->parent)) {
      y = z->parent->parent->right;
      if (is_red(y)) {
	z->parent->colour = y->colour = node::black;
	z->parent->parent->colour = node::red;
	z = z->parent->parent;
      }
      else {
	if (is_right_child(z)) { z = z->parent; rotate_left(z); }
	z->parent->colour = node::black;
	z->parent->parent->colour = node::red;
	rotate_right(z->parent->parent);
      }
    }
    else {
      y = z->parent->parent->left;
      if (is_red(y)) {
	z->parent->colour = y->colour = node::black;
	z->parent->parent->colour = node::red;
	z = z->parent->parent;
      }
      else {
	if (is_left_child(z)) { z = z->parent; rotate_right(z); }
	z->parent->colour = node::black;
	z->parent->parent->colour = node::red;
	rotate_left(z->parent->parent);
      }
    }

  root->colour = node::black;

  return z;
}

interval_tree_base::node*
interval_tree_base::first_intersecting(node* x, const interval& i) {
  if (x != &nil)  x = minimum_perhaps_intersecting(x, i);

  while (! (x == &nil || overlaps(x->first, i)))
    x = successor_perhaps_intersecting(x, i);

  return x;
}

interval_tree_base::node*
interval_tree_base::next_intersecting(node* x, const interval& i) {
  do {
    x = successor_perhaps_intersecting(x, i);
  } while (! (x == &nil || overlaps(x->first, i)));

  return x;
}

void interval_tree_base::dump(int level, char side,
			      const node* p, const node* pparent) const {
  if (p == &nil)  return;
  std::clog << p;
  for (int i = -2; i < level; i++)  std::clog << " ";
  std::clog << side << " " << (is_red(p)? "red  " : "black")
	    << " " << p->first << " max:" << p->max_zlimit+1;
  if (p->parent != pparent)  std::clog << " borked parent ptr";
  std::clog << "\n";
  dump(level + 1, 'L', p->left, p);
  dump(level + 1, 'R', p->right, p);
}

void interval_tree_base::dump_intersecting_r(node* x, const interval& i) {
  if (i.zstart() < x->left->max_zlimit)
    dump_intersecting_r(x->left, i);

  if (overlaps(x->first, i))  std::cout << ' ' << x->first;
  else  std::cout << " [" << x->first << ']';

  if (x->right != &nil && x->first.zstart() < i.zlimit())
    dump_intersecting_r(x->right, i);
}

void interval_tree_base::dump_intersecting_i(node* x, const interval& i) {
  bool go_left = true;

  while (x != &nil) {
    if (go_left)
      while (i.zstart() < x->left->max_zlimit)  x = x->left;

    if (overlaps(x->first, i))  std::cout << ' ' << x->first;
    else  std::cout << " [" << x->first << ']';

    if (x->right != &nil && x->first.zstart() < i.zlimit()) {
      x = x->right;
      go_left = true;
    }
    else {
      while (x->parent != &nil && is_right_child(x))  x = x->parent;
      x = x->parent;
      go_left = false;
    }
  }
}

} // namespace sam
