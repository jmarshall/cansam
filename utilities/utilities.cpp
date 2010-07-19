#include "utilities/utilities.h"

#include <ostream>

#include "sam/version.h"

void print_version(std::ostream& stream, const char* name) {
  stream << name << " (Cansam) " << sam::version() << "\n"
"Copyright (C) 2010 Genome Research Ltd.\n"
"This is free software: you are free to change and redistribute it.\n"
"There is NO WARRANTY, to the extent permitted by law.\n";
}
