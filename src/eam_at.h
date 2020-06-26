#ifndef EAM_AT_H
#define EAM_AT_H

/** Finnis-Sinclair form as used by Ackland and Thetford to include
 * corections for small interatomic separations.  Phil. Mag. A, 56
 * (1987) No. 1, p. 15-30.
*/

#include "object.h"
#include "potential.h"
#include "eam.h"

void eam_at_parms(POTENTIAL* object, EAM_PARMS* parms);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
