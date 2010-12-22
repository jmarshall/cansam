/*  utilities.cpp -- Various library support functions.

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

#include "lib/utilities.h"

#include <string>
#include <iomanip>

using std::string;

namespace sam {

extern const char format::hexadecimal_digits[] = "0123456789ABCDEF";

// Removes a trailing line terminator, whether it be LF, CR, or CR-LF.
// (Usually CR would be because there was a CR-LF terminator and the LF has
// already been elided.)
string& chomp(string& s) {
  string::size_type len = s.length();
  if (len >= 2 && s.compare(len - 2, 2, "\x0d\x0a") == 0)
    s.erase(len - 2);
  else if (len >= 1 && (s[len - 1] == '\x0d' || s[len - 1] == '\x0a'))
    s.erase(len - 1);

  return s;
}

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

make_string& make_string::operator<< (unsigned char c) {
  if (isgraph(c))  buffer << c;
  else if (c == '\0')  buffer << "\\0";
  else  buffer << "\\x" << std::hex << std::setfill('0') << std::setw(2)
	       << unsigned(c);

  return *this;
}

} // namespace sam
