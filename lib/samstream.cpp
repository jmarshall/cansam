#include "sam/stream.h"

#include <cctype>

#include "sam/exception.h"
#include "sam/rawfilebuf.h"

using std::string;

namespace sam {

// As setstate(), but returns whether setstate() would have thrown an exception.
bool samstream_base::setstate_wouldthrow(iostate state) {
  try { setstate(state); }
  catch (...) { return true; }

  return false;
}


isamstream::isamstream(const std::string& filename, openmode mode) {
  cansam::rawfilebuf sbuf;
  sbuf.open(filename.c_str(), mode);
}

isamstream::isamstream(std::streambuf* sbuf, openmode mode)
  : samstream_base(sbuf, sambamio::new_in(sbuf)) {
  // FIXME read the heaaders, probably
}

isamstream::~isamstream() {
  // FIXME...
}

isamstream& isamstream::operator>> (alignment& aln) {
  try {
    if (! io->get(*this, aln))
      if (setstate_wouldthrow(eofbit | failbit) && (exceptions() & eofbit))
	throw sam::failure("eof");
  }
  catch (sam::failure& e) {
    if (setstate_wouldthrow(failbit)) {
      e.set_filename("thefilename"); // FIXME
      throw;
    }
  }
  catch (sam::exception& e) {
    if (setstate_wouldthrow(badbit)) {
      e.set_filename("thefilename"); // FIXME
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
