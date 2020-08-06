/// @file cansam/interval.h
/// Classes for intervals and sequence intervals

/*  Copyright (C) 2011-2012 Genome Research Ltd.

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

#ifndef CANSAM_INTERVAL_H
#define CANSAM_INTERVAL_H

#include <iosfwd>
#include <string>

#include <stdint.h>

#include "cansam/types.h"

namespace sam {

class alignment;

/** @class sam::interval cansam/interval.h
    @brief Interval within an unspecified sequence
*/
class interval {
public:
  /// Construct an empty interval
  interval() : zstart_(0), zlimit_(0) { }

  /// Construct a zero-based, half-open interval
  interval(coord_t zstart, coord_t end) : zstart_(zstart), zlimit_(end) { }

  /// Construct an interval from a "START-END"-style string
  explicit interval(const std::string& text) { assign(text); }

  /// Construct a copy of an interval
  interval(const interval& i) : zstart_(i.zstart_), zlimit_(i.zlimit_) { }

  //  Destroy this interval object (not interesting enough to warrant ///)
  ~interval() { }

  /// Copy an interval
  interval& operator= (const interval& i)
    { zstart_ = i.zstart_; zlimit_ = i.zlimit_; return *this; }

  /// Assign to this interval from a "START-END"-style string
  interval& assign(const std::string& text, size_t pos = 0);

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

/** @class sam::seqinterval cansam/interval.h
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

  /// Construct a seqinterval from a "NAME:START-END"-style string
  explicit seqinterval(const std::string& text) { assign(text); }

  /// Construct a seqinterval representing ALN's span on its reference
  explicit seqinterval(const alignment& aln);

  /// Construct a copy of a seqinterval
  seqinterval(const seqinterval& i) : interval(i), name_(i.name_) { }

  //  Destroy this seqinterval object (not interesting enough to warrant ///)
  ~seqinterval() { }

  /// Copy a seqinterval
  seqinterval& operator= (const seqinterval& i)
    { name_ = i.name_; zstart_ = i.zstart_; zlimit_ = i.zlimit_; return *this; }

  /// Assign to this seqinterval
  seqinterval& assign(const std::string& name, coord_t zstart, coord_t end)
    { name_ = name; zstart_ = zstart; zlimit_ = end; return *this; }

  /// Assign to this seqinterval
  seqinterval& assign(const char* name, coord_t zstart, coord_t end)
    { name_ = name; zstart_ = zstart; zlimit_ = end; return *this; }

  /// Assign to this seqinterval from a "NAME:START-END"-style string
  seqinterval& assign(const std::string& text, size_t pos = 0);

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

} // namespace sam

#endif
