#include "fixedVelocity.h"
#include "group.h"
#include "ddcMalloc.h"
#include "object.h"


typedef struct fixVelocityParms_st { int set; THREE_VECTOR v;} FIXVELOCITYPARMS; 
void fixedVelocity_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in)
{
}
void fixedVelocity_velocityUpdate(int mode, int k,GROUP *g, STATE *state, double time, double dt)
{
   FIXVELOCITYPARMS *parms = (FIXVELOCITYPARMS *)g->parm; 
   if (parms->set)
   {
      state->vx[k] = parms->v.x; 
      state->vy[k] = parms->v.y; 
      state->vz[k] = parms->v.z; 
   }
}
void fixedVelocity_parms(GROUP *gp)
{
   gp->itype = FIXEDVELOCITY;
   FIXVELOCITYPARMS *parms = (FIXVELOCITYPARMS *)ddcMalloc(sizeof(FIXVELOCITYPARMS));
   gp->parm = (void*)parms; 
   gp->write_dynamics = NULL; 
   gp->velocityUpdate= (void (*)(int,int,GROUP*,void*,double,double))fixedVelocity_velocityUpdate; 
   gp->Update= (void (*)(GROUP *, int, void *, double, double))fixedVelocity_Update; 
   parms->set = object_testforkeyword((OBJECT *)gp, "v");
   object_get((OBJECT*)gp, "v", &parms->v, WITH_UNITS,     3, "0.0 0.0 0.0",     "l/t", NULL);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
