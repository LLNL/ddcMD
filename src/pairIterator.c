#include "pairIterator.h"
#include <assert.h>
#include "ddcMalloc.h"

/** iter.ii == iter.jj is the signal that this is an end iterator.
 * Since no atom can be a pair with itself this can't happen
 * naturally. */

void pairIter_advance(PAIR_ITERATOR* iter)
{
   assert(iter->ii != iter->jj);
   iter->advance(iter);
}




/* Local Variables: */
/* tab-width: 3 */
/* End: */
