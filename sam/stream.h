/** @file  sam/stream.h
    @brief Classes for SAM/BAM input/output streams

FIXME: These aren't really streams, so it would be good to come up with
a different name.
*/

#ifndef CANSAM_STREAM_H
#define CANSAM_STREAM_H

#include <string>
#include <vector>

#include "sam/types.h"
#include "sam/alignment.h"
#include "sam/header.h"

namespace sam {

/** Flags used when opening a stream.  */
enum openmode {
  in,
  out,
  trunc,
  binary,	///< Binary format, i.e., BAM
  compressed,	///< Stream contents are (to be) compressed
  uncompressed 	///< No compression

};

class samstream_base {
protected:
  // @cond private
  struct reference {
    std::string name;
    coord_t length;
  };
  std::vector<reference> reftable;
  std::vector<header> headers;

  void read_reftable();

  class sambamio {
  public:
    static sambamio* new_in(std::streambuf* sbuf);
    static sambamio* new_out(std::streambuf* sbuf, openmode);

    virtual ~sambamio() { }

    virtual bool get_headers(samstream_base&) = 0;
    virtual bool get(samstream_base&, alignment&) = 0;
  };

  class bamio;
  class samio;
  class gzsamio;

  sambamio* io;
  // @endcond
};

/** @class sam::isamstream sam/stream.h
    @brief SAM/BAM input stream

FIXME This isn't a stream in the sense that you can stream @e anything from it;
you only get to read %sam %alignment records.  Maybe we can come up with a
better name than (by parallel with ifstream, istringstream) a "stream backed
by a sam".  */
class isamstream : public samstream_base {
public:
  isamstream(const std::string& filename, openmode mode = in);
  isamstream(std::streambuf* sbuf, openmode mode = in);
  ~isamstream();

  /// Seek back to the first %alignment record in the stream.
  isamstream& rewind();

  /// Read an %alignment.
  isamstream& get(alignment& aln) { io->get(*this, aln); return *this; }
};

/// Read an %alignment from a SAM/BAM stream
isamstream& operator>> (isamstream& stream, alignment& aln)
  { return stream.get(aln); }

/** @class sam::osamstream sam/stream.h
    @brief SAM/BAM output stream */
class osamstream : public samstream_base {
public:
  osamstream(const std::string& filename, openmode mode = out);
  osamstream(std::streambuf* sbuf, openmode mode = out);
  ~osamstream();
};

/// Write the %alignment to a SAM/BAM stream
osamstream& operator<< (osamstream& stream, const alignment& aln);

/** @brief Returns the mode flags indicated by the filename extension.
@details Returns an appropriate openmode for .bam, .%sam, and .sam.gz.
Filenames are matched case-insensitively, and unrecognised extensions are
equivalent to .%sam.  */
openmode extension(const std::string& filename);

} // namespace sam

#endif
