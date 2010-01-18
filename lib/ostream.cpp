#include <ostream>

#include "sam/header.h"
#include "sam/alignment.h"

#include "lib/utilities.h"  // FIXME NUKE-ME want make_string for dump_on()

using std::string;

namespace sam {

std::ostream& operator<< (std::ostream& out, const header& header) {
#ifdef WE_HAVE_COMMENTS_BUGGER_IT
  return out << header.str();
#else
  out << '@' << header.type();
  for (header::const_iterator it = header.begin(); it != header.end(); ++it) {
    out << '\t';
    if (it->tag()[0])
      out << it->tag() << ':';
    out << it->value_str();
  }

  return out;
#endif
}

std::ostream& operator<< (std::ostream& out, const header::tagfield& field) {
  const char* limit = header::tagfield::nexttab(field.tag_);
  return out << string(field.tag_, limit - field.tag_);
}

std::ostream& operator<< (std::ostream& out, header::const_iterator it) {
  return out << static_cast<const void*>(it.ptr);
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

  for (alignment::const_iterator it = aln.begin(); it != aln.end(); ++it)
    out << '\t' << *it;

  out.flags(flags);
  return out;
}

std::ostream& operator<< (std::ostream& out, const alignment::tagfield& aux) {
  static string bam_only = "cCsSI";

  char type = aux.type();
  if (bam_only.find(type) != string::npos)  type = 'i';
  return out << aux.tag() << ':' << type << ':' << aux.value();
}

std::ostream& operator<< (std::ostream& out, alignment::const_iterator it) {
  return out << static_cast<const void*>(it.ptr);
}

// FIXME nuke me or otherwise hide me!
void alignment::dump_on(std::ostream& out, const_iterator marker) const {
  make_string text;
  const char* s = p->data();
  const char* limit = &s[p->capacity() - sizeof(block_header)];
  while (s < limit)
    if (s == marker.ptr)  text << " [" << *s++ << "] ";
    else  text << *s++;

  out << "Capacity:" << p->capacity() << ", cindex:" << p->h.cindex
      << ", data:" << string(text) << "\n";
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
