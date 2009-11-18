/** @file   sam/alignment-old.h
    @brief  Class for SAM/BAM alignment records

    Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
*/

#ifndef CANSAM_ALIGNMENT_OLD_H
#define CANSAM_ALIGNMENT_OLD_H

#include <iosfwd>
#include <string>

#include "sam/types.h"
#include "sam/tagfield.h"

namespace sam {

#if 0
// FIXME do something about namespaces so you don't have to say sam::PAIRED etc

/** The @e flags field of an alignment contains a number of flag bits pertaining
    to the read that has been aligned and a few pertaining to its mate.  */
enum alignment_flag {
  PAIRED              = 0x001, ///< The read is paired with a mate
  PROPER_PAIRED       = 0x002, ///< The read and its mate are mapped as a proper pair
  UNMAPPED            = 0x004, ///< The read is unmapped
  MATE_UNMAPPED       = 0x008, ///< The mate is unmapped
  REVERSE_STRAND      = 0x010, ///< The read's strand (0: forward; 1: reverse)
  MATE_REVERSE_STRAND = 0x020, ///< The mate's strand (0: forward; 1: reverse)
  FIRST_IN_PAIR       = 0x040, ///< The read is the first of a pair
  SECOND_IN_PAIR      = 0x080, ///< The read is the second of a pair
  NONPRIMARY          = 0x100, ///< This %alignment is not the primary one for the read
  QUALITY_FAILED      = 0x200, ///< Failed platform or vendor quality checks
  DUPLICATE           = 0x400, ///< PCR duplicate or optical duplicate
};
#endif

/** @class sam::alignment sam/alignment.h
    @brief SAM/BAM %alignment record */
class alignment {
public:
  /// @name Constructors
  //@{
  alignment() { }  ///< Default constructor
  /// Constructor with all (non-tagged) fields
  alignment(const std::string& qname, int flags, const std::string& rname,
    coord_t pos, int mapq, const std::string& cigar,
    const std::string& mate_rname, coord_t mate_pos, scoord_t isize,
    const std::string& seq, const std::string& qual)
    : qname_(qname), flags_(flags), rname_(rname), pos_(pos), mapq_(mapq),
      cigar_(cigar), mate_rname_(mate_rname), mate_pos_(mate_pos),
      isize_(isize), seq_(seq), qual_(qual) { }
  //@}

  /// @name Field Accessors
  //@{
  std::string qname() const; ///< Query name
  std::string& qname();

  int flags() const { return flags_; } ///< Bitwise combination of alignment_flag flags
  int& flags() { return flags_; }

  std::string rname() const; ///< Reference name

  coord_t  pos() const;  ///< Leftmost position (1-based)
  coord_t& pos();        ///< Leftmost position (1-based)
  coord_t  zpos() const; ///< Leftmost position (0-based)
  coord_t& zpos();       ///< Leftmost position (0-based)

  int mapq() const { return mapq_; } ///< Mapping quality

  std::string cigar() const;

  std::string mate_rname() const; ///< Mate's reference name (for paired reads)

  coord_t  mate_pos() const;  ///< Mate's leftmost position (1-based)
  coord_t& mate_pos();        ///< Mate's leftmost position (1-based)
  coord_t  mate_zpos() const; ///< Mate's leftmost position (0-based)
  coord_t& mate_zpos();       ///< Mate's leftmost position (0-based)

  scoord_t isize() const { return isize_; }
  scoord_t& isize() { return isize_; }

  std::string seq() const;

  std::string qual() const;

  const tagfield::array& auxen() const;
  tagfield::array& auxen(); ///< Array of optional auxiliary fields
  //@}

  /// Determines whether the %alignment has a particular optional field.
  /// @param tag  The two-character field tag to be checked.
  /// @return     The field's type character ('i', 'Z', etc), or '\\0' if the
  ///             %alignment does not have a field with the given tag.
  char has(const char* tag) const;

  std::string  aux(const char* tag) const;
  std::string  aux(const char* tag, const std::string& default_value) const;
  std::string& aux(const char* tag);

  int  aux_int(const char* tag) const;
  int  aux_int(const char* tag, int default_value) const;
  int& aux_int(const char* tag);

  float  aux_float(const char* tag) const;
  float  aux_float(const char* tag, float default_value) const;
  float& aux_float(const char* tag);

  double  aux_double(const char* tag) const;
  double  aux_double(const char* tag, double default_value) const;
  double& aux_double(const char* tag);

  /// @name Derived Information
  //@{
  /// @brief The @e strand flag, as '+' or '-'.
  /// @details Returns the strand information encoded in the flags field.
  /// @return  Either '+' or '-'.
  char strand() const { return (flags_ & REVERSE_STRAND)? '-' : '+'; }

  /// @brief The @e mate-strand flag, similarly.
  /// @details Returns the mate read's strand information encoded in the
  /// flags field.  This is only meaningful if the %alignment in fact has
  /// a mate, i.e., if it is paired.
  /// @return  Either '+' or '-'.
  char mate_strand() const { return (flags_ & MATE_REVERSE_STRAND)? '-' : '+'; }

#if 0
  /// Returns the pairing order information encoded in the flags field.
  /// @return  -1 when the @e first flag is set, +1 when the @e second flag
  /// is set, or 0 if neither (or both) is set.
  int order() const
    { return order_value[(flags_ & (FIRST_IN_PAIR | SECOND_IN_PAIR)) >> 6]; }
#endif

  coord_t right_pos() const;  ///< Rightmost position (1-based)
  coord_t right_zpos() const; ///< Rightmost position (0-based)

  //@}

  static int to_flags(const char* begin, const char* end);
  static int to_flags(const std::string&);

private:
  std::string qname_;
  int flags_;
  std::string rname_;
  coord_t pos_;
  int mapq_;
  std::string cigar_;
  std::string mate_rname_;
  coord_t mate_pos_;
  scoord_t isize_;
  std::string seq_;
  std::string qual_;

  tagfield::array auxen_;

  // FIXME Setting HIDE_UNDOC_MEMBERS=YES will suppress this from the doco
  // @cond private
  friend std::istream& operator>> (std::istream&, alignment&);
  // @endcond
};

/// @brief Compares alignments by genomic location
/// @details Blah blah blah.
bool operator< (const alignment& a, const alignment& b);
inline bool operator> (const alignment& a, const alignment& b) { return b < a; }

std::istream& operator>> (std::istream&, alignment&);

/** Writes an %alignment record to a %stream as text in SAM format,
    @b without a trailing newline character.  */
std::ostream& operator<< (std::ostream& strm, const alignment& aln);

} // namespace sam

#ifndef DONT_USE_THEM
using sam::PAIRED;
using sam::PROPER_PAIRED;
using sam::UNMAPPED;
using sam::MATE_UNMAPPED;
// etc
//    sam::PAIRED, sam::PROPER_PAIRED, sam::UNMAPPED, sam::MATE_UNMAPPED,
//    sam::REVERSE_STRAND, sam::MATE_REVERSE_STRAND, sam::FIRST_IN_PAIR, 
//    sam::SECOND_IN_PAIR, sam::NONPRIMARY, sam::QUALITY_FAILED, sam::DUPLICATE
#endif

#endif
