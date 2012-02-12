/*  utilities.cpp -- Support routines common to the various utilities.

    Copyright (C) 2010-2011 Genome Research Ltd.

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

#ifndef UTILITIES_UTILITIES_H
#define UTILITIES_UTILITIES_H

#include <iosfwd>
#include <string>

// Prints Cansam version number as PROGNAME's version number and
// brief copyright and (lack of) warranty information to STREAM.
void print_version(std::ostream& stream, const char* progname);

// Returns whether standard input appears not to have been redirected,
// i.e., is coming from the user.  This lets us determine when we should
// read copiously from standard input (when it's redirected from a BAM file)
// or when it would be more useful to produce e.g. a usage display (when the
// user probably doesn't want to type a SAM file by hand); for example, when
// a utility is invoked with no arguments.
bool cin_likely_from_user();

// Returns PATH with any leading directories and trailing extensions removed.
std::string basename(const std::string& path);

#endif
