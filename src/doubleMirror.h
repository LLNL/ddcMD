#ifndef DOUBLE_MIRROR_H
#define DOUBLE_MIRROR_H
#include "group.h"

void doubleMirror_Update(GROUP* g, int mode, void* state, double time_in, double dt_in);
void doubleMirror_velocityUpdate(int mode, int k, GROUP* g, STATE *state, double time, double dt);
void doubleMirror_velocityUpdate1(int mode, int k, GROUP* g, SPECIES *species, THREE_VECTOR V, THREE_VECTOR F, THREE_VECTOR R, double time, double dt);
void doubleMirror_parms(GROUP* gp);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
