/** @file  sam/stream.h
    @brief Classes for SAM/BAM input/output streams

FIXME: These aren't really streams, so it would be good to come up with
a different name.
*/

#ifndef CANSAM_STREAM_H
#define CANSAM_STREAM_H

#include <ios>
#include <string>
#include <vector>

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

};
#endif

class samstream_base : public std::ios {
public:
  // FIXME But wait!  Who's responsible for the io?
  virtual ~samstream_base() { }

  /// The filename associated with this stream, or empty if none or unknown
  std::string filename() const { return filename_; }

  /// Set an associated filename
  void set_filename(const std::string& filename) { filename_ = filename; }

protected:
  class sambamio;
  samstream_base(std::streambuf* sbuf, samstream_base::sambamio* io0)
    : std::ios(sbuf), io(io0) { }

  samstream_base() { }

  // @cond private
  bool setstate_wouldthrow(iostate state);

  class sambamio {
  public:
    static sambamio* new_in(std::streambuf* sbuf);
    static sambamio* new_out(std::streambuf* sbuf, openmode);

    virtual ~sambamio() { }

    // FIXME maybe these should be taking the collection in fact
    virtual bool get_headers(samstream_base&) = 0;

    // Returns true when an alignment has been successfully read,
    // false at EOF, or throws an exception on formatting or I/O errors.
    virtual bool get(samstream_base&, alignment&) = 0;

    virtual void put(samstream_base&, const alignment&) = 0;
  };

  class bamio;
  class samio;

  sambamio* io;
  // @endcond

private:
  std::string filename_;
};

/** @class sam::isamstream sam/stream.h
    @brief SAM/BAM input stream

FIXME This isn't a stream in the sense that you can stream @e anything from it;
you only get to read sam alignment records.  Maybe we can come up with a
better name than (by parallel with ifstream, istringstream) a "stream backed
by a sam".  */
class isamstream : public samstream_base {
public:
  isamstream(const std::string& filename, openmode mode = in);
  isamstream(std::streambuf* sbuf, openmode mode = in);
  ~isamstream();

  /// Read an alignment record
  /** Similarly to @a std::istream's extraction operators, this reads one
  alignment record into @a aln and returns the stream.  The @a iostate flags
  will be set at EOF or upon errors, and exceptions will be thrown accordingly
  as selected via @a ios::exception().  */
  isamstream& operator>> (alignment& aln);

  /// Seek back to the first alignment record in the stream
  isamstream& rewind();

  // FIXME  Some form of seek/tell -- or maybe that's in samstream_base
};

/** @class sam::osamstream sam/stream.h
    @brief SAM/BAM output stream */
class osamstream : public samstream_base {
public:
  osamstream(const std::string& filename, openmode mode = out);
  osamstream(std::streambuf* sbuf, openmode mode = out);
  ~osamstream();

  /// Write the alignment
  osamstream& operator<< (const alignment& aln);
};


/** @brief Returns the mode flags indicated by the filename extension.
@details Returns an appropriate openmode for .bam, .sam, and .sam.gz.
Filenames are matched case-insensitively, and unrecognised extensions are
equivalent to .sam.  */
std::ios_base::openmode extension(const std::string& filename);

} // namespace sam

#endif
