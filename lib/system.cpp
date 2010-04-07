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

#include "sam/streambuf.h"

namespace sam {

unsigned long get_vss() {
  rawfilebuf statfile("/proc/self/stat", O_RDONLY);
  if (! statfile.is_open())
    return 0;

  return info.virtual_size;
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
