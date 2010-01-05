/** @file  sam/stream.h
    @brief Classes for SAM/BAM input/output streams */

#ifndef CANSAM_STREAM_H
#define CANSAM_STREAM_H

#include <ios>
#include <string>

#include "sam/types.h"
#include "sam/alignment.h"
#include "sam/header.h"

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

/** @class sam::samstream_base sam/stream.h
    @brief Base class for SAM/BAM streams

This is a base class for the SAM/BAM stream hierarchy; it is unlikely to be
usefully instantiated itself.

@note This hierarchy doesn't provide a stream in the sense that @e anything can
be streamed to/from it; it only accepts SAM alignment records.  The name does
not parallel @c fstream and @c stringstream; it does not mean "a stream backed
by a sam", whatever that might be.  */
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
  class sambamio;

  // @cond private
  bool setstate_wouldthrow(iostate state);

  samstream_base(std::streambuf* sbuf, bool owned, sambamio* io0)
    : std::ios(sbuf), io(io0), filename_(), owned_rdbuf_(owned) { }

  sambamio* io;
  // @endcond

private:
  class bamio;
  class samio;

  std::string filename_;
  bool owned_rdbuf_;
};

/** @class sam::isamstream sam/stream.h
    @brief SAM/BAM input stream
*/
class isamstream : virtual public samstream_base {
public:
  /// Construct an input stream by opening a file
  explicit isamstream(const std::string& filename, openmode mode = in);

  /// Construct an input stream associated with an already-opened stream buffer
  explicit isamstream(std::streambuf* sbuf);

  virtual ~isamstream();

  /// Read an alignment record
  /** Similarly to @c std::istream's extraction operators, this reads one
  alignment record into @a aln and returns the stream.  The @c iostate flags
  will be set at EOF or upon errors, and exceptions will be thrown accordingly
  as selected via @c exceptions().  */
  isamstream& operator>> (alignment& aln);

  /// Seek back to the first alignment record in the stream
  isamstream& rewind();

  // FIXME  Some form of seek/tell -- or maybe that's in samstream_base
};

/** @class sam::osamstream sam/stream.h
    @brief SAM/BAM output stream */
class osamstream : virtual public samstream_base {
public:
  /// Construct an output stream by opening a file
  explicit osamstream(const std::string& filename, openmode mode = out);

  /// Construct an output stream associated with an already-opened stream buffer
  explicit osamstream(std::streambuf* sbuf);

  virtual ~osamstream();

  /// Write the alignment
  osamstream& operator<< (const alignment& aln);
};

/** @class sam::samstream sam/stream.h
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
