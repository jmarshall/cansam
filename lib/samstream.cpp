/*  samstream.cpp -- Classes for SAM/BAM input/output streams.

    Copyright (C) 2010-2012 Genome Research Ltd.

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

#include "cansam/sam/stream.h"

#include <fstream>
#include <cctype>    // for tolower()

#include <unistd.h>  // for STDIN_FILENO etc

// FIXME do we actually need all these, or could we just forward declare?
#include "cansam/sam/alignment.h"
#include "cansam/exception.h"
#include "cansam/streambuf.h"
#include "lib/sambamio.h"

using std::string;

namespace sam {

bool samstream_base::is_open() const {
  if (sam::streambuf* sbuf = dynamic_cast<sam::streambuf*>(rdbuf()))
    return sbuf->is_open();
  else if (std::filebuf* fbuf = dynamic_cast<std::filebuf*>(rdbuf()))
    return fbuf->is_open();
  else
    return rdbuf()? true : false;
}

void samstream_base::close() {
  if (sam::streambuf* sbuf = dynamic_cast<sam::streambuf*>(rdbuf()))
    sbuf->close();
  else if (std::filebuf* fbuf = dynamic_cast<std::filebuf*>(rdbuf()))
    fbuf->close();

  // FIXME maybe this should delete rdbuf() (if owned_rdbuf_) and/or io
}

samstream_base::~samstream_base() {
  delete io;

  if (owned_rdbuf_) {
    exceptions(goodbit);  // Prevent exceptions caused by emptying rdbuf().
    std::streambuf* sbuf = rdbuf(NULL);
    delete sbuf;
  }
}

// As per setstate(), but returns whether setstate() would have thrown an
// exception, i.e., whether a state listed in exceptions() has been set.
bool samstream_base::setstate_wouldthrow(iostate state) {
  try { setstate(state); }
  catch (...) { return true; }

  return false;
}

std::streambuf*
new_and_open(const std::string& filename, std::ios::openmode mode) {
  // TODO Eventually might look for URL schemes and make a different streambuf

  if (filename == "-") {
    rawfilebuf* sbuf = new rawfilebuf();
    // TODO On some platforms, may need setmode(O_BINARY) or similar
    sbuf->attach((mode & std::ios::out)? STDOUT_FILENO : STDIN_FILENO);
    return sbuf;
  }
  else
    return new rawfilebuf(filename.c_str(), mode & ~compressed);

  // FIXME "& ~compressed" is because compressed == ate... hmmm... that's a hack
}

isamstream::isamstream(const std::string& filename, openmode mode)
  : samstream_base(new_and_open(filename, mode | in), true) {
  set_filename(filename);
  if (is_open())
    io = sambamio::new_in(*this);
}

isamstream::isamstream(std::streambuf* sbuf)
  : samstream_base(sbuf, false) {
  io = sambamio::new_in(*this);
}

isamstream::~isamstream() {
  // FIXME...
}

#if 0
isamstream& isamstream::rewind() {
  // FIXME
  return *this;
}
#endif


/* These operators set iostate bits in response to exceptions from io->get()
and from the underlying streambuf, and propagate the exceptions if instructed
to do so by exceptions().

The formatting layers throw sam::bad_format exceptions, so these cause failbit
to be set; Cansam's streambufs specifically throw sam::* exceptions other than
sam::bad_format, and other streambufs might throw anything at all but surely
not sam::bad_format, so other exceptions set badbit, corresponding to serious
streambuf problems such as I/O errors.  */

isamstream& isamstream::operator>> (collection& headers) {
  try {
    if (! io) throw "hmmmm"; // FIXME
    io->get(*this, headers);
  }
  catch (const sambamio::eof_exception&) {
    // An exception was thrown while setting eofbit within io->get().
    // Propagate it as an externally-known exception type.
    throw std::ios::failure("eof");
  }

  return *this;
}

isamstream& isamstream::operator>> (alignment& aln) {
  try {
    if (! io->get(*this, aln)) {
      // There are no more records, so set failbit but don't throw (thus
      // differing from standard streams) as this is not a failure as such.
      // Leave eofbit as is, as it means EOF-seen-on-streambuf rather than
      // no-more-records, and it is almost certainly already set anyway.
      (void) setstate_wouldthrow(failbit);
    }
  }
  catch (const sambamio::eof_exception&) {
    // An exception was thrown while setting eofbit within io->get().
    // Propagate it as an externally-known exception type.
    throw std::ios::failure("eof");
  }
  catch (sam::bad_format& e) {
    // If an exception is to be thrown, annotate and rethrow the original
    // exception (rather than the generic one thrown by plain setstate()).
    if (setstate_wouldthrow(failbit)) {
      e.set_filename(filename());
      throw;
    }
  }
  catch (sam::exception& e) {
    if (setstate_wouldthrow(badbit)) {
      e.set_filename(filename());
      throw;
    }
  }
  catch (...) {
    // FIXME  Possibly worth propagating this as a sam::exception or so,
    // so that our callers need expect only sam::* exceptions.
    if (setstate_wouldthrow(badbit))
      throw;
  }

  return *this;
}


osamstream::osamstream(const std::string& filename, openmode mode)
  : samstream_base(new_and_open(filename, mode | out), true) {
  set_filename(filename);
  if (is_open())
    io = sambamio::new_out(mode);
}

osamstream::~osamstream() {
  // FIXME  What exactly should this condition be?
  // Probably close() above should flush too
  // or maybe ~samio() should do the flushing

  if (is_open() && io)
    io->flush(*this);
  // FIXME catch...
}

osamstream& osamstream::operator<< (const collection& headers) {
  io->put(*this, headers);
  // FIXME exceptions...
  return *this;
}

osamstream& osamstream::operator<< (const alignment& aln) {
  try {
    io->put(*this, aln);
  }
  // FIXME hmmm...
  catch (sam::exception& e) {
    if (setstate_wouldthrow(badbit)) {
      e.set_filename(filename());
      throw;
    }
  }
  catch (...) {
    if (setstate_wouldthrow(badbit))
      throw;
  }

  return *this;
}


std::ios_base::openmode extension(const string& filename) {
  using std::ios;

  string::size_type dotpos = filename.rfind('.');
  if (dotpos != string::npos) {
    // If the extension has fewer than 4 characters (i.e., is shorter
    // than ".bam" or ".sam"), look at the double-extension instead, as
    // we must be looking for ".sam.gz".
    if (filename.length() - dotpos < 4)
      dotpos = filename.rfind('.', dotpos - 1);
  }

  if (dotpos == string::npos)
    return sam_format;

  string ext(filename, dotpos);
  for (string::size_type i = 0; i < ext.length(); i++)
    ext[i] = ::tolower(ext[i]);

  if (ext == ".bam")  return bam_format;
  else if (ext == ".sam")  return sam_format;
  else if (ext == ".sam.gz")  return sam_format | compressed;
  else  return sam_format;
}

} // namespace sam
