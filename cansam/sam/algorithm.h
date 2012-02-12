/// @file cansam/sam/algorithm.h
/// Algorithms operating on alignment records

/*  Copyright (C) 2010-2012 Genome Research Ltd.

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

#ifndef CANSAM_SAM_ALGORITHM_H
#define CANSAM_SAM_ALGORITHM_H

#include <functional>

#include "cansam/sam/alignment.h"

namespace sam {

/// Function object class for comparing alignments
class less_by_qname :
  public std::binary_function<const alignment&, const alignment&, bool> {
public:
  /// Returns whether the query name of @a a is lexicographically less than
  /// that of @a b
  bool operator() (const alignment& a, const alignment& b) const
    { return cmp_by_qname(a, b) < 0; }
};

/// Function object class for comparing alignments
class equal_by_qname :
  public std::binary_function<const alignment&, const alignment&, bool> {
public:
  /// Returns whether the query name of @a a is the same as that of @a b
  bool operator() (const alignment& a, const alignment& b) const
    { return cmp_by_qname(a, b) == 0; }
};

/// Function object class for hashing alignments
class hash_by_qname : public std::unary_function<const alignment&, size_t> {
public:
  /// Returns a hash value derived from the alignment's query name
  size_t operator() (const alignment& aln) const {
    size_t result = 0;
    for (const char* s = aln.qname_c_str(); *s; s++)
      result = (result * 131) + *s;
    return result;
  }
};

} // namespace sam

#endif
