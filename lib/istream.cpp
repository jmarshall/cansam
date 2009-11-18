#include <istream>
#include <string>

#include "sam/header.h"
#include "sam/alignment.h"

#include "lib/utilities.h"

using std::string;

namespace sam {

int peek(std::istream& strm) {
  int c = strm.peek();
  return (c == EOF || c == '@')? c : 'A';
}

std::istream& operator>> (std::istream& strm, header& hdr) {
  if (strm.peek() != '@') {
    // It's not a header, but could be EOF, an alignment record, or badly
    // formatted data.  In all these cases, set failbit to indicate the error.
    strm.setstate(std::ios::failbit);
    return strm;
  }

  string line;
  if (! getline(strm, line))
    return strm;

  hdr.assign(chomp(line));
  return strm;
}

#if 0
// FIXME nuke me
coord_t to_int(const string& str,
	       string::size_type begin, string::size_type end) {
  const char* data = str.data();
  const char* s = data + begin;
  const char* lim = data + end;

  coord_t val = 0;
  while (s < lim)
    val = 10 * val + *s++ - '0';
  return val;
}
#endif

std::istream& operator>> (std::istream& strm, alignment& aln) {
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

  int c = strm.peek();
  if (c == EOF || c == '@') {
    strm.setstate(std::ios::failbit);
    return strm;
  }

  string line;
  if (! getline(strm, line))
    return strm;

  aln.assign(chomp(line));
  return strm;
#if 0
  string::size_type tabpos, pos;
  
  tabpos = line.find('\t');
  aln.qname_.assign(line, 0, tabpos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.flags_ = alignment::to_flags(line.substr(pos, tabpos - pos));

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.rname_.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.pos_ = to_int(line, pos, tabpos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.mapq_ = to_int(line, pos, tabpos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.cigar_.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  if (line.compare(pos, tabpos - pos, "=") == 0)
    aln.mate_rname_ = aln.rname_;
  else
    aln.mate_rname_.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.mate_pos_ = to_int(line, pos, tabpos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  if (line[pos] == '-')  aln.isize_ = -to_int(line, pos + 1, tabpos);
  else  aln.isize_ = to_int(line, pos, tabpos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.seq_.assign(line, pos, tabpos - pos);

  pos = tabpos + 1;
  tabpos = line.find('\t', pos);
  aln.qual_.assign(line, pos, tabpos - pos);

  // FIXME clear() calls a whole lot of destructors -- better to assign the
  // first min(oldcount,newcount) tagfields and then resize(newcount) or
  // push_back(ctor) any remaining new ones, rather than dtor x oldcount
  // followed by ctor x newcount.
  aln.auxen_.clear();
  while (tabpos != string::npos) {
    pos = tabpos + 1;
    tabpos = line.find('\t', pos);
    aln.auxen_.push_back(tagfield(line.substr(pos, tabpos - pos)));
  }
#endif
}

} // namespace sam
