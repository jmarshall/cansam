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
class refsequence;

/** @class sam::collection sam/collection.h
    @brief Header information for a collection of SAM/BAM records */
class collection {
private:
  // TODO  This might be becoming boost::ptr_vector<header>
  typedef std::vector<header*> header_array;

public:
  /// Construct an empty collection
  collection();

  /// Construct a copy of a collection, by copying all the headers within
  collection(const collection& collection);

  //  Destroy this collection object (not interesting enough to warrant ///)
  ~collection();

  /// Copy a collection, by copying all the headers within
  collection& operator= (const collection& collection);

  // @cond infrastructure
  typedef header_array::iterator iterator;  // FIXME or something...
  typedef header_array::const_iterator const_iterator;
  typedef header_array::reference reference;
  typedef header_array::const_reference const_reference;
  typedef header_array::difference_type difference_type;
  // @endcond

  iterator begin() { return headers.begin(); }
  const_iterator begin() const { return headers.begin(); }

  iterator end() { return headers.end(); }
  const_iterator end() const { return headers.end(); }

  void push_back(const std::string& header_line);

  bool empty() const { return headers.empty(); }
  void clear();

  refsequence& findseq(const std::string& name);
  refsequence& findseq(const char* name);
  refsequence& findseq(int index);

  // FIXME prob not public
  static collection& find(unsigned cindex) { return *collections[cindex]; }

private:
  int cindex;
  std::vector<header*> headers;
  std::vector<refsequence*> refseqs;
  std::map<std::string, int> refmap;

  // FIXME maybe just make these local to lib/collection.cpp
  static std::vector<collection*> collections;
};

} // namespace sam

#endif
