#ifndef BERENDSEN_H
#define BERENDSEN_H
#include "state.h"
#include "group.h"

void berendsen_Update(GROUP* g, int mode, void* state, double time_in, double dt_in);
void berendsen_velocityUpdate(int mode, int k, GROUP* g, STATE *state, double time, double dt);
void berendsen_parms(GROUP* gp);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
