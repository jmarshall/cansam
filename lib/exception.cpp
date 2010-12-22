/*  exception.cpp -- Exceptions thrown by Cansam classes and functions.

    Copyright (C) 2010 Genome Research Ltd.

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

#include "sam/exception.h"

#include <sstream>
#include <cstring>  // for strerror()

namespace sam {

namespace {

// Returns whether S ends in a non-space character, i.e., is a complete phrase.
bool complete_phrase(const std::string& s) {
  return (! s.empty()) && s[s.length() - 1] != ' ';
}

} // unnamed namespace

// Returns "message [[for ] "filename"] [[at record] N]"
const char* bad_format::what() const throw() {
  std::ostringstream s;

  s << sam::exception::what();
  if (! filename().empty()) {
    if (complete_phrase(s.str()))  s << " for ";
    s << '"' << filename() << '"';
  }

  if (recnum() != 0) {
    if (complete_phrase(s.str()))  s << " at record ";
    s << recnum();
  }

  what_text_ = s.str();
  return what_text_.c_str();
}

// Returns "message [[for ] "filename"]: strerror
const char* system_error::what() const throw() {
  what_text_ = sam::exception::what();

  if (! filename().empty()) {
    if (complete_phrase(what_text_))  what_text_ += " for ";
    what_text_ += '"';
    what_text_ += filename();
    what_text_ += '"';
  }

  what_text_ += ": ";
  what_text_ += strerror(errnum());

  return what_text_.c_str();
}

} // namespace sam
