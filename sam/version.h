/// @file sam/version.h
/// Cansam library version information

/*  Copyright (C) 2010 Genome Research Ltd.

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

#ifndef CANSAM_VERSION_H
#define CANSAM_VERSION_H

#include <string>

/** @file
The Cansam library version number is a 3-part @c major.minor.patch number,
provided both as a preprocessor constant and a text string.  While the former
denotes the Cansam headers your code is compiled against and the latter the
Cansam library used at link- or run-time, in most circumstances both will
represent the same value.

The preprocessor macro #CANSAM_VERSION encodes the version as an integer
literal in the same way as <a href="http://www.boost.org/">Boost</a>'s
version macro, with two decimal digits of @c patch and three of @c minor.
The version() function returns the version number in text form, with the
final @c ".nn" omitted if @c patch is 0.

This header file is not included by other Cansam header files, so your code
should @c @#include it itself where necessary.

@note While the utilities supplied with the Cansam library use this library
version number directly as their own version number, external code should not
-- as such code is not updated in concert with the library.  */

/// Library version, as an integer @c Mmmmpp
#define CANSAM_VERSION  690

namespace sam {

/// Library version, as @c "n.nnn[.nn]"
std::string version();

} // namespace sam

#endif
