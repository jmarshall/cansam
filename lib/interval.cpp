/*  interval.cpp -- Classes for intervals and sequence intervals.

    Copyright (C) 2011-2012 Genome Research Ltd.

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

#include "cansam/interval.h"

#include <limits>
#include <ostream>

#include "cansam/exception.h"
#include "lib/utilities.h"

using std::string;

namespace sam {

namespace {

int32_t parse_numeral(const char*& s0, int32_t default_value) {
  const char* s = s0;
  int32_t value = 0;

  while (true)
    if (*s >= '0' && *s <= '9')  value = 10 * value + (*s++ - '0');
    else if (*s == ',')  s++;
    else  break;

  if (s > s0)  s0 = s;
  else  value = default_value;

  return value;
}

} // unnamed namespace


interval& interval::assign(const string& text, size_t pos) {
  // The string is either "[START]", "[START]-[END]", or "[START]+[LENGTH]".
  const char* s = text.c_str() + pos;

  zstart_ = parse_numeral(s, 1) - 1;

  switch (*s) {
  case '-':
    s++;
    zlimit_ = parse_numeral(s, std::numeric_limits<int32_t>::max());
    break;

  case '+':
    s++;
    zlimit_ = zstart_ + parse_numeral(s, 0);
    break;

  case '\0':
    zlimit_ = zstart_ + 1;
    break;

  default:
    break;
  }

  if (*s != '\0')
    throw bad_format(make_string()
	<< "Invalid interval value ('" << text.substr(pos) << "')");

  return *this;
}

seqinterval& seqinterval::assign(const string& text, size_t pos) {
  size_t colonpos = text.rfind(':');
  if (colonpos < pos)  colonpos = string::npos;

  name_.assign(text, pos, colonpos - pos);
  if (colonpos == string::npos)  interval::assign("-");
  else  interval::assign(text, colonpos+1);

  return *this;
}

std::ostream& operator<< (std::ostream& stream, const interval& i) {
  return stream << i.start() << '-' << i.end();
}

std::ostream& operator<< (std::ostream& stream, const seqinterval& i) {
  return stream << i.name() << ':' << i.start() << '-' << i.end();
}

} // namespace sam
