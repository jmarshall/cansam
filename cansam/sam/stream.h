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
class exception;
class sambamio;

/** @class sam::samstream_base cansam/sam/stream.h
    @brief Base class for SAM/BAM streams

This is a base class for the SAM/BAM stream hierarchy; it cannot be
instantiated itself.

Because it is not just characters being transferred, there is buffering and
other state in the <tt>[io]samstream</tt> object as well as the associated
@c streambuf.  Hence using the inherited @c rdbuf(streambuf*) method to change
the associated stream buffer has an undefined effect.

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

  /// Set initial exceptions mask for subsequent samstream objects
  /** By default, each newly-constructed SAM/BAM stream object has an
  exceptions mask of @c failbit|badbit, so throws exceptions on all formatting
  errors and I/O errors, including failure to open files.  (This differs from
  standard iostreams, which throw no exceptions by default.)

  This function can be used to change this default for all
  subsequently-constructed stream objects.  (Use the usual @c exceptions()
  method to change an individual stream's exception mask.)  */
  static void initial_exceptions(iostate except);

protected:
  // @cond infrastructure
  void setstate_maybe_rethrow(iostate state);
  void setstate_maybe_rethrow(iostate state, sam::exception& exception);

  samstream_base();
  samstream_base(std::streambuf* sbuf, bool owned);

  void reset_closed_or_throw();
  void reset_closed_or_throw(std::streambuf* sbuf, bool owned);

  void open_into_rdbuf(const std::string& filename, openmode mode);
  virtual void close_() = 0;

  sambamio* io;
  // @endcond

private:
  friend class sambamio;

  std::string filename_;
  bool owned_rdbuf_;

  static iostate initial_exceptions_;
};

/** @class sam::isamstream cansam/sam/stream.h
    @brief SAM/BAM input stream
*/
class isamstream : virtual public samstream_base {
public:
  /// Construct an unopened input stream
  isamstream() { }

  /// Construct an input stream by opening a file
  explicit isamstream(const std::string& filename);

  /// Construct an input stream associated with an already-opened stream buffer
  explicit isamstream(std::streambuf* sbuf);

  virtual ~isamstream();

  /// Open a file
  void open(const std::string& filename);

  /// Associate with an already-opened stream buffer
  void open(std::streambuf* sbuf);

  /// Read the collection of headers
  /** Similarly to @c std::istream's extraction operators, this reads @b all
  header records into @a headers and returns the stream.  The @c iostate flags
  will be set at EOF or upon errors, and exceptions will be thrown accordingly
  as selected via @c exceptions().

  Header records appear at the start of a file, so this must be the first
  extraction operator used.

  The stream retains a reference to @a headers and hands it to
  subsequently-extracted alignment records; therefore @a headers's lifetime
  must end after the lifetimes of the stream and any resulting @c alignment
  objects.  */
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

protected:
  // @cond infrastructure
  virtual void close_();
  // @endcond
};

/** @class sam::osamstream cansam/sam/stream.h
    @brief SAM/BAM output stream */
class osamstream : virtual public samstream_base {
public:
  /// Construct an unopened output stream
  osamstream() { }

  /// Construct an output stream by opening a file
  explicit osamstream(const std::string& filename, openmode mode = out);

  /// Construct an output stream associated with an already-opened stream buffer
  explicit osamstream(std::streambuf* sbuf, openmode mode = out);

  virtual ~osamstream();

  /// Open a file
  void open(const std::string& filename, openmode mode = out);

  /// Associate with an already-opened stream buffer
  void open(std::streambuf* sbuf, openmode mode = out);

  /// Write the collection of headers
  osamstream& operator<< (const collection& headers);

  /// Write the alignment
  osamstream& operator<< (const alignment& aln);

  /// Apply @a manipulator to this stream
  osamstream& operator<< (ios_base& (*manipulator)(ios_base&))
    { manipulator(*this); return *this; }

  /// Flush any uncommitted output
  osamstream& flush();

protected:
  // @cond infrastructure
  virtual void close_();
  // @endcond
};

/** @class sam::samstream cansam/sam/stream.h
    @brief SAM/BAM input/output stream  */
class samstream : public isamstream, public osamstream {
public:
  /// Construct an unopened stream
  samstream() { }

  /// Construct a stream by opening a file
  explicit samstream(const std::string& filename, openmode mode = in|out);

  /// Construct a stream associated with an already-opened stream buffer
  explicit samstream(std::streambuf* sbuf, openmode mode = in|out);

  virtual ~samstream();

  /// Open a file
  void open(const std::string& filename, openmode mode = in|out);

  /// Associate with an already-opened stream buffer
  void open(std::streambuf* sbuf, openmode mode = in|out);

protected:
  // @cond infrastructure
  virtual void close_();
  // @endcond
};

/// Returns the mode flags indicated by the filename extension
/** Returns an appropriate openmode for .bam, .sam, and .sam.gz.
Filenames are matched case-insensitively, and unrecognised extensions are
equivalent to .sam.  */
std::ios_base::openmode extension(const std::string& filename);

} // namespace sam

#endif
