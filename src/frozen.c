#include "group.h"
#include <assert.h>
#include "ddcMalloc.h"
#include "state.h"
#include "three_algebra.h"
#include "particle.h"
#include "object.h"
#include "mpiTypes.h"

/** The purpose of the frozen group is to define a group of atoms that
 *  will not move during the simulation.  Furthermore, we want to
 *  accomplish this without disturbing the velocities of these atoms
 *  (such as by setting them to zero).  To accomplish this, we set the
 *  velocity of the atoms to zero in the FRONT_TIMESTEP part of the
 *  velocityUpdate function, but we also store a copy of the value so
 *  that we can restore it in the BACK_TIMESTEP velocityUpdate.  This
 *  way the velocity is always zero when the integrator tries to use it
 *  to move particles, but is maintained at its original value
 *  everywhere else.
 *
 *  Note that since assignment may be called between the front and back
 *  time steps it is necessary to register the velocity storage array so
 *  that the data will move with the particles.  It is true that the
 *  frozen particles don't move, but the local particle order could
 *  change as mobile particles move out of the domain and there is the
 *  possibility that the domain centers do move so even stationary
 *  particles might change domains.
 *
 */


typedef struct frozen_parms_st { THREE_VECTOR *vLast; } FROZEN_PARMS; 
void frozen_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in)
{
}

void frozen_velocityUpdate(int mode, int k, GROUP *g, STATE *state, double time, double dt)
{
   double *vx = state->vx; 
   double *vy = state->vy; 
   double *vz = state->vz; 
   FROZEN_PARMS *parm = (FROZEN_PARMS *)(g->parm); 
   THREE_VECTOR *vLast = parm->vLast; 
   switch (mode)
   {
     case FRONT_TIMESTEP:
      vLast[k].x = vx[k];  vx[k] = 0;
      vLast[k].y = vy[k];  vy[k] = 0;
      vLast[k].z = vz[k];  vz[k] = 0;
      break;
     case BACK_TIMESTEP:
      vx[k] = vLast[k].x;
      vy[k] = vLast[k].y;
      vz[k] = vLast[k].z;
      break;
     default:
      assert(0==1);
   }
}
void frozen_parms(GROUP *gp)
{
   gp->itype = FROZEN;
   gp->write_dynamics = NULL; 
   gp->velocityUpdate= (void (*)(int,int,GROUP*,void*,double,double))frozen_velocityUpdate; 
   gp->Update= (void (*)(GROUP *, int, void *, double, double))frozen_Update; 
   FROZEN_PARMS *parm= ddcMalloc(sizeof(FROZEN_PARMS)); 
   parm->vLast= NULL; 
   gp->parm = (void *)parm; 
   
   particleRegisterinfo((void*)(&parm->vLast), sizeof(*parm->vLast), sizeof(*parm->vLast), threeVector_MPIType(),NULL,NULL);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
