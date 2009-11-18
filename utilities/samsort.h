/** @file   utilities/samsort.h
    @brief  Infrastructure for extending the @c samsort utility
*/
#ifndef SAMSORT_H
#define SAMSORT_H

namespace sam { class alignment; }

/** @class  alignment_comparator utilities/samsort.h
    @brief  Helper class for adding new comparison functions to the
	    @c samsort utility.

New code can be added to @c samsort simply by adding your object file to
@c SAMSORT_OBJS in the makefile.  Your source file should include a static
object informing the existing code of your new comparison function, like this:

@code
bool less_than(const sam::alignment& a, const sam::alignment& b) {
  // ...
}

static const alignment_comparator
  my_comparator("custom", "Order by my custom criteria", less_than);
  @endcode

This object's constructor runs before @e main() to add your function to
the list; there is no need to change the existing @c samsort code itself
to inform it of your new function.  */
class alignment_comparator {
public:
  /// Signature for comparison functions, which should return true iff @a a @< @a b
  typedef bool compare(const sam::alignment& a, const sam::alignment& b);

  /** @brief Constructor, registering the comparison function.
  */
  alignment_comparator(const char* name, const char* description,
		       compare* function);

  const char* description; ///< The description that was provided
  compare* comparer; ///< The comparison function that was provided
};

#endif
