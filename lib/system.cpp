/*  system.cpp -- Miscellaneous system-specific utility functions.

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

#if defined __APPLE__

#include <mach/task.h>
#include <mach/mach_init.h>

namespace sam {

unsigned long get_vss() {
  struct task_basic_info info;
  mach_msg_type_number_t count = TASK_BASIC_INFO_COUNT;

  if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&info, &count)
      != KERN_SUCCESS)
    return 0;

  return info.virtual_size;
}

} // namespace sam

#elif defined __linux__

#include <fcntl.h>

#include "cansam/streambuf.h"

namespace sam {

unsigned long get_vss() {
  rawfilebuf statfile;
  if (! statfile.open("/proc/self/stat", O_RDONLY))
    return 0;

  return 0;  // FIXME
}

} // namespace sam

#else

#include <sys/resource.h>

namespace sam {

unsigned long get_vss() {
  struct rusage usage;

  if (getrusage(RUSAGE_SELF, &usage) != 0)
    return 0;

  return usage.ru_maxrss;
}

} // namespace sam
#endif
