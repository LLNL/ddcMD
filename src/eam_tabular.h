#ifndef EAM_TABULAR_H
#define EAM_TABULAR_H

/** Support for generalized tabular eam potentials.  If you can tabulate
 *  the pair energy, pair density and embedding functions then this form
 *  will work for you.  */

#include "object.h"
#include "potential.h"
#include "eam.h"

void eam_tabular_parms(POTENTIAL* object, EAM_PARMS* parms);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
