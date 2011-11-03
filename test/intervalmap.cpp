/*  test/intervalmap.cpp -- Tests for intervals and interval containers.

    Copyright (C) 2011 Genome Research Ltd.

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

#include <iostream>

#include "test/test.h"
#include "sam/intervalmap.h"

void search(sam::interval_multimap<char>& m, const sam::seqinterval& i) {
  m.dump_intersecting(i);

  std::cout << i << ':';
  sam::interval_multimap<char>::iterator_pair range = m.intersecting_range(i);
  for (sam::interval_multimap<char>::iterator it = range.first;
       it != range.second; ++it)
    std::cout << ' ' << it->first;
  std::cout << "  iterator!\n\n";
}

void test_intervals(test_harness& t) {
  sam::seqinterval si("X", 1000, 5000);
  std::cout << si << "\n";

  sam::interval_multimap<char> repeats;

  repeats.insert(std::make_pair(si, 'A'));

  sam::seqinterval si2("X", 4000, 8000);
  repeats.insert(std::make_pair(si2, 'B'));

  sam::seqinterval si3("X", 400, 3000);
  repeats.insert(std::make_pair(si3, 'C'));

  repeats.insert(std::make_pair(sam::seqinterval("X", 800, 6000), 'D'));
  repeats.insert(std::make_pair(sam::seqinterval("X", 3000, 4200), 'E'));

  std::clog << "\nNow the example from the book...\n\n";

  repeats.insert(std::make_pair(sam::seqinterval("Y", 16-1, 22), 'a'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 8-1, 10), 'b'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 15-1, 24), 'c'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 5-1, 9), 'd'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 25-1, 31), 'e'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 17-1, 20), 'f'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 19-1, 21), 'g'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 6-1, 11), 'h'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 26-1, 27), 'i'));
  repeats.insert(std::make_pair(sam::seqinterval("Y", 1-1, 4), 'j'));

  search(repeats, sam::seqinterval("Y", 0, 50));
  search(repeats, sam::seqinterval("Y", 8, 19));
  search(repeats, sam::seqinterval("Y", 12, 20));
  search(repeats, sam::seqinterval("Y", 19, 24));
  search(repeats, sam::seqinterval("X", 3999, 5000));

  sam::interval_multimap<char>::iterator_pair range =
    repeats.intersecting_range(sam::seqinterval("X", 1, 50000));

  for (sam::interval_multimap<char>::iterator it = range.first;
       it != range.second; ++it)
    std::clog << "X:" << it->first << " -> " << it->second << "\n";

  range = repeats.intersecting_range(sam::seqinterval("Y", 1, 50000));

  for (sam::interval_multimap<char>::iterator it = range.first;
       it != range.second; ++it)
    std::clog << "Y:" << it->first << " -> " << it->second << "\n";
}
