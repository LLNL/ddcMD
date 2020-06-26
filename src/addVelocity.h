#ifndef ADDVELOCITY_H
#define ADDVELOCITY_H

#include "transform.h"

/** The addVelocity transform is multipurpose.  When it is specified
 * with type=ADDVELOCITY it will add a specified constant to the
 * particle velocities on a per group or per species basis.  When it is
 * specified with type=SETVELOCITY it will set the center of mass
 * velocity of the system as a whole.  The function that implements the
 * transform (addVelocity) and the parms structure is the same either
 * way.  The only difference is the routine that is called to get the
 * parmeters from the object.  */

/** Set up parameters for an ADDVELOCITY transform */
void* addVelocity_parms(TRANSFORM*);
/** Set up parameters for a SETVELOCITY transform */
void* setVelocity_parms(TRANSFORM*);
/** Perform *either* an ADDVELOCITY or SETVELOCITY transform. */
void addVelocity(TRANSFORM*);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
