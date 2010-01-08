/** @file  sam/collection.h
    @brief Class for SAM/BAM header collective infos
*/

#ifndef CANSAM_COLLECTION_H
#define CANSAM_COLLECTION_H

#include <string>
#include <vector>
#include <map>

#include "sam/types.h"

namespace sam {

#if 0
/*. @class sam::reference sam/collection.h
    @brief Reference sequence record, corresponding to an @@SQ header */
class reference {
public:
  reference(const std::string& name, coord_t length)
    : name_(name), length_(length) { }

  ~reference() { }

  /// Name of the sequence
  std::string name() const { return name_; }

  /// Sequence length
  coord_t length() const { return length_; }

  void set_name(const std::string& name) { name_ = name; }
  void set_length(coord_t length) { length_ = length; }

private:
  std::string name_;
  coord_t length_;
};
#endif

/** @class sam::collection sam/collection.h
    @brief Header information for a collection of SAM/BAM records */
class collection {
public:
  collection();
  ~collection();

  // FIXME or call it rindex or so?
  int findseq(const std::string& rname) const;

  std::string rname(int rindex) const;

#if 0
  // @cond private
  struct reference {
    std::string name;
    coord_t length;
  };
  std::vector<reference> reftable;
  std::vector<header> headers;

  void read_reftable();
#endif

  // FIXME prob not public
  static collection& find(unsigned cindex) { return *collections[cindex]; }

private:
  //std::vector<something> refseqs;
  std::map<std::string, int> refmap;

  // FIXME maybe just make this local to lib/collection.cpp
  static std::vector<collection*> collections;
};

} // namespace sam

#endif
