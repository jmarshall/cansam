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

class header;

/** @class sam::collection sam/collection.h
    @brief Header information for a collection of SAM/BAM records */
class collection {
public:
  collection();
  collection(const collection& collection);
  collection& operator= (const collection& collection);
  ~collection();

  // @cond infrastructure
  typedef std::vector<header*>::iterator iterator;  // FIXME or something...

  // @endcond

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
  std::vector<header*> headers;
  //std::vector<something> refseqs;
  std::map<std::string, int> refmap;

  // FIXME maybe just make this local to lib/collection.cpp
  static std::vector<collection*> collections;
};

} // namespace sam

#endif
