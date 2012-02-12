/*  collection.cpp -- Class for a set of SAM/BAM headers.

    Copyright (C) 2010-2012 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
 3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
    nor the names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND ITS CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */

#include "cansam/sam/header.h"
#include "cansam/exception.h"
#include "lib/sambamio.h"  // for push_back() flags

#include "lib/utilities.h"

using std::string;

namespace sam {

collection::collection() {
  allocate_cindex();
  refseqs_in_headers = false;
}

collection::~collection() {
  clear();
  free_cindex();
}

std::vector<collection*>
  collection::collections(1, static_cast<collection*>(NULL));

void collection::allocate_cindex() {
  // TODO Reuse old ones after a while
  cindex = collections.size();
  collections.push_back(this);
}

template <typename InputIterator>
void delete_each(InputIterator first, InputIterator last) {
  for (InputIterator it = first; it != last; ++it)
    delete *it;
}

void collection::clear() {
  refnames.clear();
  rgroups.clear();

  if (! refseqs_in_headers)
    delete_each(refseqs.begin(), refseqs.end());
  refseqs.clear();

  delete_each(headers.begin(), headers.end());
  headers.clear();

  refseqs_in_headers = false;
}

// FIXME Make me class static?
static refsequence unmapped_refseq("*", 0, -1);

refsequence& collection::findseq_(int index) const {
  if (index >= 0 && index < int(refseqs.size()))
    return *refseqs[index];
  else if (index == -1)
    return unmapped_refseq;
  else
    throw sam::exception("Reference sequence index out of range"); // FIXME
}

refsequence& collection::findseq_(const string& name) const {
  refname_map::const_iterator it = refnames.find(name);
  if (it != refnames.end())
    return *(it->second);
  else if (name == "*")  // FIXME Is "*" going to be in refnames, or not?
    return unmapped_refseq;
  else
    throw sam::exception(make_string()
	<< "No such reference sequence ('" << name << "')");
}

refsequence& collection::findseq_(const char* name) const {
  if (name[0] == '*' && name[1] == '\0')
    return unmapped_refseq;

  refname_map::const_iterator it = refnames.find(string(name));
  if (it != refnames.end())
    return *(it->second);
  else
    throw sam::exception(make_string()
	<< "No such reference sequence ('" << name << "')");
}

readgroup& collection::findgroup_(const std::string& id) const {
  readgroup_map::const_iterator it = rgroups.find(id);
  if (it == rgroups.end())
    throw sam::exception(make_string()
	<< "No such read group ('" << id << "')");

  return *(it->second);
}

void collection::push_back(const std::string& header_line) {
  std::string text = header_line;
  size_t pos = 0;
  while ((pos = text.find('\t', pos)) != string::npos)
    text[pos++] = '\0';

  // FIXME  Is all-flags-always correct?
  push_back(text, add_header | add_refseq | add_refname);
}

// TEXT is NUL-delimited.
void collection::push_back(const std::string& text, int flags) {
  header* hdr;
  if (text.compare(0, 3, "@SQ") == 0) {
    refsequence* rhdr =
	new refsequence(text, (flags & add_refseq)? refseqs.size() : -1);

    if (flags & add_refseq)   refseqs.push_back(rhdr);
    if (flags & add_refname)  refnames[rhdr->name()] = rhdr;
    hdr = rhdr;
  }
  else if (text.compare(0, 3, "@RG") == 0) {
    readgroup* rghdr = new readgroup(text);
    rgroups[rghdr->id()] = rghdr;
    hdr = rghdr;
  }
  else
    hdr = new header(text);

  if (flags & add_header)  headers.push_back(hdr);
}

} // namespace sam
