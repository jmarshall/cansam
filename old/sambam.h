#ifndef SAMBAM_H
#define SAMBAM_H

#include <iosfwd>
#include <istream>
#include <string>

using std::string;

typedef long coord_t;
typedef long scoord_t;

enum SamFlag {
  PAIRED              = 0x001,
  PROPER_PAIRED       = 0x002,
  UNMAPPED            = 0x004,
  MATE_UNMAPPED       = 0x008,
  REVERSE_STRAND      = 0x010,
  MATE_REVERSE_STRAND = 0x020,
  FIRST_IN_PAIR       = 0x040,
  SECOND_IN_PAIR      = 0x080,
  NONPRIMARY          = 0x100,
  QUALITY_FAILED      = 0x200,
  DUPLICATE           = 0x400,
  };

class SamAlignment {
public:
  SamAlignment() { }
  SamAlignment(const string& qname0, int flags0, const string& rname0,
    coord_t pos0, int mapq0, const string& cigar0, const string& mate_rname0,
    coord_t mate_pos0, scoord_t isize0, const string& seq0, const string& qual0)
    : qname(qname0), flags(flags0), rname(rname0), pos(pos0), mapq(mapq0),
      cigar(cigar0), mate_rname(mate_rname0), mate_pos(mate_pos0),
      isize(isize0), seq(seq0), qual(qual0) { }

  string qname;
  int flags;
  string rname;
  coord_t pos;
  int mapq;
  string cigar;
  string mate_rname;
  coord_t mate_pos;
  scoord_t isize;
  string seq;
  string qual;
  // Note extras starts with a tab if it is non-empty
  string extras; // FIXME Chop optional extra fields up and implement properly

  char strand() const { return (flags & REVERSE_STRAND)? '-' : '+'; }
  char mate_strand() const { return (flags & MATE_REVERSE_STRAND)? '-' : '+'; }
};

std::istream& operator>> (std::istream&, SamAlignment&);
std::ostream& operator<< (std::ostream&, const SamAlignment&);

class SamAlignmentFromFile {
public:
  string qname() const { return string(delim[0] + 1, delim[1] - delim[0]); }
  int flags() const;
  string rname() const;
  coord_t pos() const;
  int mapq() const;
  string cigar() const;
  string mate_rname() const;
  coord_t mate_pos() const;
  scoord_t isize() const;
  string seq() const;
  string qual() const;

  bool is_paired() const { return flags() & PAIRED; }
  bool is_proper_paired() const { return flags() & PROPER_PAIRED; }
  bool is_unmapped() const { return flags() & UNMAPPED; }
  bool mate_is_unmapped() const { return flags() & MATE_UNMAPPED; }
  bool is_first() const { return flags() & FIRST_IN_PAIR; }
  bool is_second() const { return flags() & SECOND_IN_PAIR; }
  bool is_secondary() const { return flags() & NONPRIMARY; }
  bool is_duplicate() const { return flags() & DUPLICATE; }

  char strand() const { return (flags() & REVERSE_STRAND)? '-' : '+'; }
  char mate_strand() const {return (flags() & MATE_REVERSE_STRAND)? '-' : '+';}

  friend std::istream& operator>> (std::istream&, SamAlignmentFromFile&);
  friend std::ostream& operator<< (std::ostream&, const SamAlignmentFromFile&);

private:
  friend class SamInput;

  string line;
  const char* delim[30];
};

class SamInput {
public:
  bool getAlignment(SamAlignment& aln);
  bool getAlignment(SamAlignmentFromFile& aln);

private:
  std::istream is;
};

#endif

/* vi:set sw=2: */
