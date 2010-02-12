#include <ostream>

#include "sam/alignment.h"
#include "sam/header.h"
#include "lib/utilities.h"  // FIXME NUKE-ME want make_string for dump_on()

using std::string;

namespace sam {

std::ostream& operator<< (std::ostream& out, const header& header) {
#if 1
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

std::ostream& operator<< (std::ostream& out, const collection& headers) {
  for (collection::const_iterator it = headers.begin();
       it != headers.end(); ++it)
    out << **it << '\n';

  if (out.flags() & std::ios::showpoint) {
    int i = 0;
    out << "Reflist:";
    for (std::vector<refsequence*>::const_iterator it = headers.refseqs.begin();
	 it != headers.refseqs.end(); ++it, ++i)
      out << " " << i << "->" << (*it)->name();
    out << "\nRefmap:";
    for (std::map<std::string, refsequence*>::const_iterator
	 it = headers.refnames.begin(); it != headers.refnames.end(); ++it)
      out << " " << it->first << "->" << it->second->name();
    out << "\n";
  }

  return out;
}

std::ostream& operator<< (std::ostream& out, const header::tagfield& field) {
#ifdef STILL_TAB_DELIMITED
  const char* limit = header::tagfield::next(field.tag_);
  return out << string(field.tag_, limit - field.tag_);
#else
  return out << field.tag_;
#endif
}

std::ostream& operator<< (std::ostream& out, header::const_iterator it) {
  return out << static_cast<const void*>(it.ptr);
}

std::ostream& operator<< (std::ostream& out, const alignment& aln) {
  std::ios_base::fmtflags flags = out.flags();

  out << aln.qname() << '\t' << std::showbase << aln.flags()
      << '\t' << aln.rname() << '\t' << std::dec << aln.pos()
      << '\t' << aln.mapq() << '\t' << aln.cigar();

  if (aln.mate_rindex() == aln.rindex() && aln.rindex() >= 0)  out << "\t=";
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
  const char* limit = &s[p->capacity()];
  while (s < limit) {
    if (s == qname_c_str())  text << "]NAME:[";
    if (s == p->data() + p->cigar_offset())  text << "]CIG:[";
    if (s == seq_raw_data())  text << "]SEQ:[";
    if (s == qual_raw_data())  text << "]QUAL:[";
    if (s == begin().ptr)  text << "]AUXEN:[";

    if (s == marker.ptr)  text << " [" << *s++ << "] ";
    else  text << *s++;
  }

  out << "Capacity:" << p->capacity() << ", cindex:" << p->h.cindex
      << ", data:[" << string(text) << "]\n";
}

} // namespace sam
