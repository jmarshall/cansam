/*  samstream.cpp -- Classes for SAM/BAM input/output streams.

    Copyright (C) 2010-2012, 2018 Genome Research Ltd.

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
#include <errno.h>

#include <unistd.h>  // for STDIN_FILENO etc

// FIXME do we actually need all these, or could we just forward declare?
#include "cansam/sam/alignment.h"
#include "cansam/exception.h"
#include "cansam/streambuf.h"
#include "lib/sambamio.h"

using std::string;

namespace sam {

namespace {

// Used as a placeholder rdbuf() for closed streams.  This class's methods are
// never invoked, as closed_io always accompanies closed_buf and doesn't use
// its rdbuf().  (We can't use NULL as a placeholder as that causes badbit to
// be set and likely throws an exception inconveniently.)
class nullbuf : public std::streambuf {
public:
  nullbuf() { }
  virtual ~nullbuf() { }
} closed_buf;

// Used as a placeholder sambamio for closed streams, throwing an exception if
// any of the stream's insertion/extraction operators are used -- rather than
// crashing as would happen if we used a NULL sambamio for closed streams.
class nullio : public sambamio {
public:
  nullio(const char* what) : error(what) { }
  virtual ~nullio() { }

  virtual bool get(isamstream&, collection&) { throw error; }
  virtual bool get(isamstream&, alignment&)  { throw error; }
  virtual void put(osamstream&, const collection&) { throw error; }
  virtual void put(osamstream&, const alignment&)  { throw error; }
  virtual void flush(osamstream&) { throw error; }

protected:
  virtual size_t xsgetn(isamstream&, char*, size_t) { throw error; }

private:
  std::logic_error error;
} closed_io("samstream is not open");

} // anonymous namespace


// For use while handling an exception.  Sets error state as per setstate(),
// but if an exception is to be thrown (due to setting a state listed in
// exceptions()), rethrows the original exception instead of std::ios::failure.
void samstream_base::setstate_maybe_rethrow(iostate state) {
  bool rethrow = false;

  try { setstate(state); }
  catch (...) { rethrow = true; }

  if (rethrow)  throw;
}

// For use while handling an exception, which must be provided as an argument.
// As above, also annotating the rethrown exception with the stream's details.
void samstream_base::setstate_maybe_rethrow(iostate state, sam::exception& e) {
  bool rethrow = false;

  try { setstate(state); }
  catch (...) { rethrow = true; }

  if (rethrow) {
    e.set_filename(filename());
    throw;
  }
}


std::ios::iostate samstream_base::initial_exceptions_ = failbit | badbit;

void samstream_base::initial_exceptions(iostate except) {
  initial_exceptions_ = except;
}


// The usual base constructor sets up a tidily-closed stream, so that if
// further construction of an open stream fails, our destructor has a known
// state from which to determine whether to deallocate  io  and  rdbuf().
// (Whenever  rdbuf() is &closed_buf,  owned_rdbuf_  will be false.)
samstream_base::samstream_base()
  : std::ios(&closed_buf), io(&closed_io), filename_(), owned_rdbuf_(false) {
  exceptions(initial_exceptions_);
}

// When the stream buffer is already known, we can short-circuit that
// and construct our  std::ios  base with the final stream buffer.
samstream_base::samstream_base(std::streambuf* sbuf, bool owned)
  : std::ios(sbuf), io(&closed_io), filename_(), owned_rdbuf_(owned) {
  exceptions(initial_exceptions_);
}

// These reset the closed stream to the same state as the corresponding
// constructors, or throw if it is already open.  (A closed stream is
// either tidily-closed or has a  non-closed_buf rdbuf()  left behind by
// a failed  open().  These functions tidy up the latter before resetting.)
void samstream_base::reset_closed_or_throw() {
  reset_closed_or_throw(&closed_buf, false);
}

void samstream_base::reset_closed_or_throw(std::streambuf* sbuf, bool owned) {
  if (io != &closed_io)  throw std::logic_error("samstream is already open");

  std::streambuf* prevbuf = rdbuf(sbuf);
  if (owned_rdbuf_)  delete prevbuf;
  owned_rdbuf_ = owned;

  filename_.clear();
}


bool samstream_base::is_open() const {
  return io != &closed_io;
}

void samstream_base::close()
try {
  if (io == &closed_io)  throw std::logic_error("samstream is already closed");

  // Invoke class-specific actions, e.g., flushing osamstreams.
  // (Derived classes' close_() methods assume their stream is properly open,
  // so all callers must check  io  or  is_open()  before invoking them.)
  close_();

  // If the stream buffer can be closed, do so.  We close it explicitly
  // rather than just via its destructor so that exceptions are propagated,
  // which is usually the point of calling this close() explicitly.
  if (sam::streambuf* sbuf = dynamic_cast<sam::streambuf*>(rdbuf()))
    sbuf->close();
  else if (std::filebuf* fbuf = dynamic_cast<std::filebuf*>(rdbuf()))
    fbuf->close();

  // Return this stream to an unopened state.
  delete io;
  io = &closed_io;

  std::streambuf* prevbuf = rdbuf(&closed_buf);
  if (owned_rdbuf_)  delete prevbuf;
  owned_rdbuf_ = false;

  filename_.clear();
}
catch (sam::bad_format& e) { setstate_maybe_rethrow(failbit, e); }
catch (sam::exception& e)  { setstate_maybe_rethrow(badbit, e); }
catch (...) { setstate_maybe_rethrow(badbit); }

samstream_base::~samstream_base() {
  // Derived classes' destructors will already have invoked their
  // class-specific close_() actions as required.

  if (io != &closed_io)  delete io;
  if (owned_rdbuf_)  delete rdbuf();
}

// Opens the specified stream into an owned rdbuf() and sets filename as
// appropriate.  If a streambuf was allocated, then it will be in rdbuf()
// and owned_rdbuf_ will be true, even if an exception is thrown thereafter;
// hence samstream_base's destructor will be able to deallocate it reliably.
void samstream_base::
open_into_rdbuf(const std::string& filename, openmode mode) {
  // TODO Eventually might look for URL schemes and make a different streambuf

  if (filename == "-") {
    if ((mode & in) && (mode & (out|app)))
      throw sam::exception("can't open standard input/output for update");

    filename_ = (mode & in)? "standard input" : "standard output";

    rawfilebuf* sbuf = new rawfilebuf;
    rdbuf(sbuf);
    owned_rdbuf_ = true;

    // TODO On some platforms, may need setmode(O_BINARY) or similar
    sbuf->attach((mode & in)? STDIN_FILENO : STDOUT_FILENO);
  }
  else {
    filename_ = filename;

    rawfilebuf* sbuf = new rawfilebuf;
    rdbuf(sbuf);
    owned_rdbuf_ = true;

    if (! sbuf->open(filename.c_str(), mode & ~compressed))
      throw sam::system_error((mode & in)? "can't open " : "can't write to ",
			      filename, errno);
  }

  // FIXME "& ~compressed" is because compressed == ate... hmmm... that's a hack
}


// Input streams
// =============

isamstream::isamstream(const std::string& filename) {
  try {
    open_into_rdbuf(filename, in | binary);
    io = sambamio::new_in(*this);
  }
  catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
  catch (...) { setstate_maybe_rethrow(failbit); }
}

isamstream::isamstream(std::streambuf* sbuf) : samstream_base(sbuf, false) {
  try {
    io = sambamio::new_in(*this);
  }
  catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
  catch (...) { setstate_maybe_rethrow(failbit); }
}

void isamstream::open(const std::string& filename)
try {
  reset_closed_or_throw();
  open_into_rdbuf(filename, in | binary);
  io = sambamio::new_in(*this);
}
catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
catch (...) { setstate_maybe_rethrow(failbit); }

void isamstream::open(std::streambuf* sbuf)
try {
  reset_closed_or_throw(sbuf, false);
  io = sambamio::new_in(*this);
}
catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
catch (...) { setstate_maybe_rethrow(failbit); }

void isamstream::close_() {
}

isamstream::~isamstream() {
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

If so-instructed, setting eofbit within io->get() produces an eof_exception
which is propagated as an externally-known exception type, std::ios::failure.

The formatting layers throw sam::bad_format exceptions, so these cause failbit
to be set; Cansam's streambufs specifically throw sam::* exceptions other than
sam::bad_format, and other streambufs might throw anything at all but surely
not sam::bad_format, so other exceptions set badbit, corresponding to serious
streambuf problems such as I/O errors.

When an exception is to be thrown, the original exception is annotated (if
it is derived from sam::exception) and rethrown rather than allowing plain
setstate() to throw a generic exception.  */

isamstream& isamstream::operator>> (collection& headers) {
  try {
    io->get(*this, headers);
  }
  catch (sambamio::eof_exception&) { throw failure("eof"); }
  catch (sam::bad_format& e) { setstate_maybe_rethrow(failbit, e); }
  catch (sam::exception& e)  { setstate_maybe_rethrow(badbit, e); }
  catch (...) { setstate_maybe_rethrow(badbit); }

  return *this;
}

isamstream& isamstream::operator>> (alignment& aln) {
  try {
    if (! io->get(*this, aln)) {
      // There are no more records, so set failbit but don't throw (thus
      // differing from standard streams) as this is not a failure as such.
      // Leave eofbit as is, as it means EOF-seen-on-streambuf rather than
      // no-more-records, and it is almost certainly already set anyway.
      try { setstate(failbit); }
      catch (...) { }
    }
  }
  catch (sambamio::eof_exception&) { throw failure("eof"); }
  catch (sam::bad_format& e) { setstate_maybe_rethrow(failbit, e); }
  catch (sam::exception& e)  { setstate_maybe_rethrow(badbit, e); }
  catch (...) { setstate_maybe_rethrow(badbit); }

  return *this;
}


// Output streams
// ==============

osamstream::osamstream(const std::string& filename, openmode mode) {
  try {
    open_into_rdbuf(filename, mode | out);
    io = sambamio::new_out(mode);
  }
  catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
  catch (...) { setstate_maybe_rethrow(failbit); }
}

osamstream::osamstream(std::streambuf* sbuf, openmode mode)
  : samstream_base(sbuf, false) {
  try {
    io = sambamio::new_out(mode);
  }
  catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
  catch (...) { setstate_maybe_rethrow(failbit); }
}

void osamstream::open(const std::string& filename, openmode mode)
try {
  reset_closed_or_throw();
  open_into_rdbuf(filename, mode | out);
  io = sambamio::new_out(mode);
}
catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
catch (...) { setstate_maybe_rethrow(failbit); }

void osamstream::open(std::streambuf* sbuf, openmode mode)
try {
  reset_closed_or_throw(sbuf, false);
  io = sambamio::new_out(mode);
}
catch (sam::exception& e)  { setstate_maybe_rethrow(failbit, e); }
catch (...) { setstate_maybe_rethrow(failbit); }

void osamstream::close_() {
  io->flush(*this);
}

osamstream::~osamstream() {
  try { if (is_open())  close_(); }
  catch (...) { }
}

/* Similarly to isamstream's operators, these methods set iostate bits in
response to exceptions from io->put() and from the underlying streambuf,
and propagate the exceptions if instructed to do so by exceptions().

This is an output stream, so sam::bad_format is unlikely but nonetheless
translated to failbit.  All other exceptions indicate serious streambuf
problems such as I/O errors, so cause badbit to be set.  */

osamstream& osamstream::operator<< (const collection& headers) {
  try {
    io->put(*this, headers);
  }
  catch (sam::bad_format& e) { setstate_maybe_rethrow(failbit, e); }
  catch (sam::exception& e)  { setstate_maybe_rethrow(badbit, e); }
  catch (...) { setstate_maybe_rethrow(badbit); }

  return *this;
}

osamstream& osamstream::operator<< (const alignment& aln) {
  try {
    io->put(*this, aln);
  }
  catch (sam::bad_format& e) { setstate_maybe_rethrow(failbit, e); }
  catch (sam::exception& e)  { setstate_maybe_rethrow(badbit, e); }
  catch (...) { setstate_maybe_rethrow(badbit); }

  return *this;
}

osamstream& osamstream::flush() {
  try {
    io->flush(*this);
  }
  catch (sam::bad_format& e) { setstate_maybe_rethrow(failbit, e); }
  catch (sam::exception& e)  { setstate_maybe_rethrow(badbit, e); }
  catch (...) { setstate_maybe_rethrow(badbit); }

  return *this;
}


// Infrastructure
// ==============

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
