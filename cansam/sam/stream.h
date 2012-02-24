/// @file cansam/sam/stream.h
/// Classes for SAM/BAM input/output streams

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

#ifndef CANSAM_SAM_STREAM_H
#define CANSAM_SAM_STREAM_H

#include <ios>
#include <string>

namespace sam {

/*. @name Additional openmode flags */
//.{
/// Consider the stream to be compressed
const std::ios_base::openmode compressed = std::ios_base::ate;

/// Flags appropriate for a SAM file
const std::ios_base::openmode sam_format
  = std::ios_base::binary ^ std::ios_base::binary; // i.e., zero

/// Flags appropriate for a BAM file (equivalent to @c binary|compressed)
const std::ios_base::openmode bam_format = std::ios_base::binary | compressed;
//.}

class alignment;
class collection;
class sambamio;

/** @class sam::samstream_base cansam/sam/stream.h
    @brief Base class for SAM/BAM streams

This is a base class for the SAM/BAM stream hierarchy; it is unlikely to be
usefully instantiated itself.

@note This hierarchy doesn't provide a stream in the sense that @e anything
can be streamed to/from it; it only accepts SAM headers and alignment records.
The name does not parallel @c fstream and @c stringstream; it does not mean
"a stream backed by a sam", whatever that might be.  */
class samstream_base : public std::ios {
public:
  virtual ~samstream_base();

  /// Returns whether the stream has been successfully opened
  bool is_open() const;

  /// Close the stream (if it is open)
  void close();

  /// The filename associated with this stream, or empty if none or unknown
  std::string filename() const { return filename_; }

  /// Set an associated filename
  void set_filename(const std::string& filename) { filename_ = filename; }

protected:
  // @cond infrastructure
  bool setstate_wouldthrow(iostate state);

  samstream_base(std::streambuf* sbuf, bool owned)
    : std::ios(sbuf), io(), filename_(), owned_rdbuf_(owned) { }

  sambamio* io;
  // @endcond

private:
  friend class sambamio;

  std::string filename_;
  bool owned_rdbuf_;
};

/** @class sam::isamstream cansam/sam/stream.h
    @brief SAM/BAM input stream
*/
class isamstream : virtual public samstream_base {
public:
  /// Construct an input stream by opening a file
  explicit isamstream(const std::string& filename, openmode mode = in);

  /// Construct an input stream associated with an already-opened stream buffer
  explicit isamstream(std::streambuf* sbuf);

  virtual ~isamstream();

  /// Read the collection of headers
  /** Blah blah blah FIXME */
  isamstream& operator>> (collection& headers);

  /// Read an alignment record
  /** Similarly to @c std::istream's extraction operators, this reads one
  alignment record into @a aln and returns the stream.  The @c iostate flags
  will be set at EOF or upon errors, and exceptions will be thrown accordingly
  as selected via @c exceptions().  */
  isamstream& operator>> (alignment& aln);

#if 0
  /// Seek back to the first alignment record in the stream
  // FIXME Or to the start of the stream, i.e., the collection?
  isamstream& rewind();
#endif

  // FIXME  Some form of seek/tell -- or maybe that's in samstream_base
};

/** @class sam::osamstream cansam/sam/stream.h
    @brief SAM/BAM output stream */
class osamstream : virtual public samstream_base {
public:
  /// Construct an output stream by opening a file
  explicit osamstream(const std::string& filename, openmode mode = out);

  /// Construct an output stream associated with an already-opened stream buffer
  explicit osamstream(std::streambuf* sbuf);

  virtual ~osamstream();

  /// Write the collection of headers
  osamstream& operator<< (const collection& headers);

  /// Write the alignment
  osamstream& operator<< (const alignment& aln);

  /// Apply @a manipulator to this stream
  osamstream& operator<< (ios_base& (*manipulator)(ios_base&))
    { manipulator(*this); return *this; }
};

/** @class sam::samstream cansam/sam/stream.h
    @brief SAM/BAM input/output stream  */
class samstream : public isamstream, public osamstream {
public:
  /// Construct a stream by opening a file
  explicit samstream(const std::string& filename, openmode mode = in|out);

  /// Construct a stream associated with an already-opened stream buffer
  explicit samstream(std::streambuf* sbuf);

  virtual ~samstream();
};

/// Returns the mode flags indicated by the filename extension
/** Returns an appropriate openmode for .bam, .sam, and .sam.gz.
Filenames are matched case-insensitively, and unrecognised extensions are
equivalent to .sam.  */
std::ios_base::openmode extension(const std::string& filename);

} // namespace sam

#endif