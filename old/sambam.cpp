#include "sambam.h"

int to_flags(const string& str) {
  int val = 0;

  string::size_type len = str.length();
  const char* s = str.data();

  if (len == 0) {
    // The field is empty, which is mostly invalid and unlikely,
    // so returning zero is not much worse than anything else.
  }
  else if (s[0] == '0' && len >= 2 && (s[1] == 'x' || s[1] == 'X')) {
    for (string::size_type i = 2; i < len; i++, s++) {
      val *= 16;
      if (*s >= 'a')  val += s[i] - 'a';
      else if (*s >= 'A')  val += s[i] - 'A';
      else  val += s[i] - '0';
    }
  }
  else if (s[0] >= '0' && s[1] <= '9') {
    int base = (s[0] == '0')? 8 : 10;
    for (string::size_type i = 0; i < len; i++, s++) {
      val *= base;
      val += *s - '0';
    }
  }
  else {
    // FIXME One day there'll be a defined text representation
  }

  return val;
}

coord_t to_int(const string& str) {
  coord_t val = 0;

  string::size_type len = str.length();
  const char* s = str.data();
  for (string::size_type i = 0; i < len; i++, s++) {
    val *= 10;
    val += s[i] - '0';
  }

  return val;
}

bool SamInput::getAlignment(SamAlignment& aln) {
  string line;
  if (! getline(is, line))
    return false;

  if (line[0] == '@')
    return false;  // FIXME ungetline!!

  string::size_type tabpos, pos;
  
  tabpos = line.find('\t');
  aln.qname.assign(line, 0, tabpos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.flags = to_flags(line.substr(pos, tabpos - pos));

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.rname.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.pos = to_int(line.substr(pos, tabpos - pos));

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.mapq = to_int(line.substr(pos, tabpos - pos));

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.cigar.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.mate_rname.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.mate_pos = to_int(line.substr(pos, tabpos - pos));

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  if (line[pos] != '-')  aln.isize = to_int(line.substr(pos, tabpos - pos));
  else  aln.isize = -to_int(line.substr(pos + 1, tabpos - pos - 1));

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.seq.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.qual.assign(line, pos, tabpos - pos);

  if (tabpos != string::npos)
    aln.extras.assign(line, tabpos, string::npos);
  else
    aln.extras.clear();

  return true;
}

std::istream& operator>> (std::istream& strm, SamAlignment& aln) {
  // FIXME Rewrite the superduper way
#if 0
  std::istream::sentry sentry(strm, true);
  if (sentry) {
    std::streambuf* sbuf = strm.rdbuf();

    int c = sbuf->sgetc();
    if (c != EOF && c != '@') {
    }
  }
#endif

  return strm;
}

std::ostream& operator<< (std::ostream& out, const SamAlignment& aln) {
  std::ios_base::fmtflags flags = out.flags();

  out << aln.qname
      << '\t' << aln.flags << '\t' << aln.rname << '\t' << aln.mapq
      << '\t' << aln.cigar;

  if (aln.mate_rname == aln.rname)  out << "\t=";
  else  out << '\t' << aln.mate_rname;

  out << '\t' << aln.mate_pos << '\t' << aln.isize << '\t' << aln.seq
      << '\t' << aln.qual << aln.extras;

  out.flags(flags);
  return out;
}

#if 0
bool SamInput::getAlignment(SamAlignment& aln) {
  return getline(is, aln.line);
}

std::istream& operator>> (std::istream& in, SamAlignment& aln) {
  return in;
}

std::ostream& operator<< (std::ostream& out, const SamAlignment& aln) {
  std::ios_base::fmtflags flags = out.flags();

  out << aln.qname() << '\t' << aln.flags() << '\t' << aln.rname() << '\t'
      << aln.mapq() << '\t' << aln.cigar() << '\t';

  if (aln.mate_rname() == aln.rname())  out << "=\t";
  else  out << aln.mate_rname() << '\t';

  out << aln.mate_pos() << '\t' << aln.isize() << '\t' << aln.seq() << '\t'
      << aln.qual();

  // FIXME also optional fields

  return out;
}
#endif

/* vi:set sw=2: */
