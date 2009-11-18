#include <ostream>

#include "sam/header.h"
#include "sam/alignment.h"

namespace sam {

std::ostream& operator<< (std::ostream& out, const header& hdr) {
  out << '@' << hdr.type();
  for (header::const_iterator it = hdr.begin(); it != hdr.end(); ++it) {
    out << '\t';
    if (! it->tag.empty())
      out << it->tag << ':';
    out << it->value;
  }

  return out;
}

std::ostream& operator<< (std::ostream& out, const alignment& aln) {
  std::ios_base::fmtflags flags = out.flags();

  out << aln.qname() << '\t' << std::showbase << aln.flags()
      << '\t' << aln.rname() << '\t' << std::dec << aln.pos()
      << '\t' << aln.mapq() << '\t' << aln.cigar();

  if (aln.mate_rname() == aln.rname())  out << "\t=";
  else  out << '\t' << aln.mate_rname();

  out << '\t' << aln.mate_pos() << '\t' << aln.isize() << '\t' << aln.seq()
      << '\t' << aln.qual();

  // FIXME
#if 0
  for (tagfield::array::const_iterator it = aln.auxen().begin();
       it != aln.auxen().end(); ++it)
    out << '\t' << *it;
#endif

  out.flags(flags);
  return out;
}

#if 0
std::ostream& operator<< (std::ostream& out, const tagfield& field) {
  out << field.tag << ':';
  if (field.type != '@')
    out << field.type << ':';

  switch (field.type) {
  case '@':
  case 'A':
  case 'H':
  case 'Z':
    out << field.str;
    break;

  case 'c':
  case 'C':
  case 's':
  case 'S':
  case 'i':
  case 'I':
    out << field.i;
    break;

  case 'f':
    out << field.f;
    break;

  case 'd':
    out << field.d;
    break;
  }

  return out;
}
#endif

} // namespace sam
