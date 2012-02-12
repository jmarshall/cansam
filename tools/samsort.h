/// @file tools/samsort.h
/// Infrastructure for extending the @c samsort utility

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

#ifndef SAMSORT_H
#define SAMSORT_H

namespace sam { class alignment; }

/** @class  alignment_comparator tools/samsort.h
    @brief  Helper class for adding new comparison functions to the
	    @c samsort utility.

New code can be added to @c samsort simply by adding your object file to
@c SAMSORT_OBJS in the makefile.  Your source file should include a static
object informing the existing code of your new comparison function, like this:

@code
bool less_than(const sam::alignment& a, const sam::alignment& b) {
  // ...
}

static const alignment_comparator
  my_comparator("custom", "Order by my custom criteria", less_than);
  @endcode

This object's constructor runs before @e main() to add your function to
the list; there is no need to change the existing @c samsort code itself
to inform it of your new function.  */
class alignment_comparator {
public:
  /// Signature for comparison functions, which should return true iff @a a @< @a b
  typedef bool compare(const sam::alignment& a, const sam::alignment& b);

  /** @brief Constructor, registering the comparison function.
  */
  alignment_comparator(const char* name, const char* description,
		       compare* function);

  const char* description; ///< The description that was provided
  compare* comparer; ///< The comparison function that was provided
};

#endif
