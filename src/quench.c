#include "group.h"
#include "state.h"


void quench_Update(GROUP *g, int mode, STATE *state, double time, double dt)
{
}

void quench_velocityUpdate(int mode, int k, GROUP *g, STATE *state, double time, double dt)
{
   double *vx = state->vx; 
   double *vy = state->vy; 
   double *vz = state->vz; 
   double *fx = state->fx; 
   double *fy = state->fy; 
   double *fz = state->fz; 
   SPECIES **species = state->species; 
   double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;
   if (vx[k] * fx[k] < 0) vx[k] = 0;
   if (vy[k] * fy[k] < 0) vy[k] = 0;
   if (vz[k] * fz[k] < 0) vz[k] = 0;
   
   double a = dt/mass; 
   vx[k] += a*fx[k];
   vy[k] += a*fy[k];
   vz[k] += a*fz[k];
}

void quench_parms(GROUP *gp)
{
   gp->itype = QUENCH;
   gp->parm = NULL;
   gp->write_dynamics = NULL; 
   gp->velocityUpdate= (void (*)(int,int,GROUP*,void*,double,double))quench_velocityUpdate; 
   gp->Update= (void (*)(GROUP *, int, void *, double, double))quench_Update; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
