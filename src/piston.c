#include "group.h"
#include "state.h"
#include "object.h"
#include "ddcMalloc.h"
#include "eq.h"

/** This is a group to implement a piston drive.  It sets the velocity
 *  of the particles in the group to vx = vy = 0, vz = f(t), where f(t)
 *  is given in the input deck.
 *
 *  This group is very close in spirit to both EXTFORCE and
 *  FIXED_VELOCITY.  However, neither of those implmentations does
 *  exactly what we want here.  I don't think either of those groups is
 *  in use so I could change and/or generalized them to be a piston, but
 *  writing this group took only 10 minutes.  We can always combine them
 *  in the future if we determine that is a good idea.  */



typedef struct piston_parms_st
{
   EQTARGET* vzeq;
} PISTON_PARMS;



void piston_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in)
{
}
void piston_velocityUpdate(int mode, int k,GROUP *gp, STATE *state, double time, double dt)
{
   PISTON_PARMS* p = gp->parm;       
   double vz_t = (p->vzeq->function)((time+dt),(void *)p->vzeq);

   state->vx[k] = 0.0;
   state->vy[k] = 0.0;
   state->vz[k] = vz_t; 
}
void piston_parms(GROUP *gp)
{
   PISTON_PARMS *parm;
   char *string; 
   parm = ddcCalloc(1, sizeof(PISTON_PARMS));
   gp->itype = PISTON;
   object_get((OBJECT *) gp, "vz", &string, STRING, 1, NULL);
   parm->vzeq = eq_parse(string,"l/t","t");
   gp->parm = parm; 
   gp->write_dynamics = NULL;
   gp->velocityUpdate= (void (*)(int,int,GROUP*,void*,double,double))piston_velocityUpdate; 
   gp->Update= (void (*) (GROUP*, int, void *, double, double))piston_Update; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
