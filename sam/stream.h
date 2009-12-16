/** @file  sam/stream.h
    @brief Classes for SAM/BAM input/output streams

FIXME: These aren't really streams, so it would be good to come up with
a different name.
*/

#ifndef CANSAM_STREAM_H
#define CANSAM_STREAM_H

#include <ios>
#include <string>

#include "sam/types.h"
#include "sam/alignment.h"
#include "sam/header.h"

namespace sam {

#if 0
/** Flags used when opening a stream.  */
enum openmode {
  in,
  out,
  trunc,
  binary,	///< Binary format, i.e., BAM
  compressed,	///< Stream contents are (to be) compressed
  uncompressed 	///< No compression

  // Extra openmode flags wanted: compressed | uncompressed
};
#endif

/** @class sam::samstream_base sam/stream.h
    @brief Base class for SAM/BAM streams  */
class samstream_base : public std::ios {
public:
  class sambamio;

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
  // @cond private
  bool setstate_wouldthrow(iostate state);

//  friend class sambamio;

  samstream_base() : std::ios(NULL) { }
  samstream_base(std::streambuf* sbuf, sambamio* io0)
    : std::ios(sbuf), io(io0), delete_rdbuf(false) { }

  sambamio* io;
  // @endcond

private:
  std::string filename_;
  bool delete_rdbuf;
};

/** @class sam::isamstream sam/stream.h
    @brief SAM/BAM input stream

FIXME This isn't a stream in the sense that you can stream @e anything from it;
you only get to read sam alignment records.  Maybe we can come up with a
better name than (by parallel with ifstream, istringstream) a "stream backed
by a sam".  */
class isamstream : virtual public samstream_base {
public:
  explicit isamstream(const std::string& filename, openmode mode = in);
  explicit isamstream(std::streambuf* sbuf, openmode mode = in);
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
  explicit osamstream(const std::string& filename, openmode mode = out);
  explicit osamstream(std::streambuf* sbuf, openmode mode = out);
  virtual ~osamstream();

  /// Write the alignment
  osamstream& operator<< (const alignment& aln);
};

#if 0
/** @class sam::samstream sam/stream.h
    @brief SAM/BAM input/output stream  */
class samstream : public isamstream, public osamstream {
public:
  explicit samstream(const std::string& filename, openmode mode = in|out);
  explicit samstream(std::streambuf* sbuf, openmode mode = in|out);
  virtual ~samstream();
};
#endif

/** @brief Returns the mode flags indicated by the filename extension.
@details Returns an appropriate openmode for .bam, .sam, and .sam.gz.
Filenames are matched case-insensitively, and unrecognised extensions are
equivalent to .sam.  */
std::ios_base::openmode extension(const std::string& filename);

} // namespace sam

#endif
