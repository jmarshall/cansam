#include "sam/stream.h"

#include <fstream>
#include <cctype>

#include "sam/exception.h"
#include "sam/streambuf.h"
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

  // FIXME maybe this should delete rdbuf() (if delete_rdbuf) and/or io
}

samstream_base::~samstream_base() {
  delete io;

  if (delete_rdbuf) {
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


/* This constructor wants to:
    init ios(NULL)
    io = NULL
    delete_rdbuf = true
    filename_ = filename
    create an appropriate sam::streambuf   [can throw]
    rdbuf(sbuf);
    open sbuf as appropriate   [can possibly throw]
    if (is_open)
      io = new_in(...)   [can throw]
*/

isamstream::isamstream(const std::string& filename, openmode mode)
  : samstream_base(new rawfilebuf, NULL) {
  //sbuf.open(filename.c_str(), mode);
}

isamstream::isamstream(std::streambuf* sbuf, openmode mode)
  : samstream_base(sbuf, sambamio::new_in(sbuf)) {
  // FIXME read the heaaders, probably
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


#if 0
std::streamsize isamstream::rdbuf_sgetn(char* buffer, std::streamsize size) {
  std::stream n;
  try {
    n = rdbuf()->sgetn(buffer, size);
  }
  catch (sam::exception& e) {
    e.set_filename(filename());
    if (setstate_wouldthrow(badbit))  throw;
    else  throw nonpropagated();
  }
  catch (...) {
    if (setstate_wouldthrow(badbit))  throw;
    else  throw nonpropagated();
  }
}
#endif

/* These operators set iostate bits in response to exceptions from io->get()
and from the underlying streambuf, and propagate the exceptions if instructed
to do so by exceptions().

The formatting layers throw sam::failure exceptions, so these cause failbit
to be set; Cansam's streambufs specifically throw sam::* exceptions other than
sam::failure, and other streambufs might throw anything at all but presumably
not sam::failure, so other exceptions set badbit, corresponding to serious
streambuf problems such as I/O errors.  */

isamstream& isamstream::operator>> (alignment& aln) {
  try {
    if (! io->get(*this, aln)) {
      // We're at EOF, so set eofbit and failbit; but (and this differs from
      // standard streams) throw only if exceptions are requested for eofbit.
      if (setstate_wouldthrow(eofbit | failbit) && (exceptions() & eofbit))
	throw sam::failure("eof");
    }
  }
  catch (sam::failure& e) {
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

// ********************************************************************

#if 0
// FIXME
int get4() { return 37; }

// Read the BAM reference table (n_ref, l_name/name/l_ref...).
void samstream_base::read_reftable() {
  int n = get4();

  char* data; // FIXME

  reftable.resize(n);
  for (int i = 0; i < n; i++) {
    int len = get4();
    reftable[i].name.assign(data, len - 1);
    data += len;
    reftable[i].length = get4();
  }
}
#endif

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
    return ios::in; // FIXME openmode(0);

  string ext(filename, dotpos);
  for (string::size_type i = 0; i < ext.length(); i++)
    ext[i] = ::tolower(ext[i]);

  // FIXME
#if 0
  if (ext == ".bam")  return ios::binary; // FIXME binary|compressed
  else if (ext == ".sam")  return openmode(2);  // FIXME fix us too
  else if (ext == ".sam.gz")  return compressed;
  else  return openmode(0);
#endif
  return ios::in;
}

} // namespace sam
