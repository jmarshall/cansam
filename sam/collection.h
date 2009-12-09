/** @file  sam/collection.h
    @brief Class for SAM/BAM header collective infos
*/

#ifndef CANSAM_COLLECTION_H
#define CANSAM_COLLECTION_H

#include <string>
#include <vector>
#include <map>

namespace sam {

/** @class sam::collection sam/collection.h
    @brief Header information for a collection of SAM/BAM records */
class collection {
public:
  collection();
  ~collection();

  // FIXME or call it rindex or so?
  int findseq(const std::string& rname) const;

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

private:
  //std::vector<something> refseqs;
  std::map<std::string, int> refmap;

  // FIXME maybe just make this local to lib/collection.cpp
  static std::vector<collection*> collections;
};

} // namespace sam

#endif
