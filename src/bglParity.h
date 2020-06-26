#ifndef BGL_PARITY_H
#define BGL_PARITY_H

#include "simulate.h"
#include "ddc.h"

/** Functions related to memory parity checking on BG/L.  These are the
 * application specific functions.  The generic support for parity error
 * handling is in parityHandler.c */

/**  This function looks for a parity failure on *all* tasks.  To avoid
 *   deadlock it is essential that slave_parityFailure be called on all
 *   slave tasks at the same time the master tasks call parityFailure.
 */
int parityFailure(SIMULATE* simulate, DDC* ddc);
void slave_parityFailure(char* data);


#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
