#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "collection.h"
#include "ddcMalloc.h"
#include "eq.h"

typedef struct extforce_parms_st
{
  EQTARGET *Fzeq;
} EXTFORCE_PARMS;




void extforce_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in)
{
        
}
void extforce_velocityUpdate(int mode, int k,GROUP *gp, STATE *state, double time, double dt)
{
   double *vx = state->vx; 
   double *vy = state->vy; 
   double *vz = state->vz; 
   double *fx = state->fx; 
   double *fy = state->fy; 
   double *fz = state->fz; 
   SPECIES **species = state->species; 
   double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;

	double a;
	EXTFORCE_PARMS *p;
	p = gp->parm;       
	double Fz_t = (p->Fzeq->function)((time+dt),(void *)p->Fzeq);
	a = dt/mass;
  	vx[k] += a*fx[k] ; 
	vy[k] += a*fy[k] ; 
  	vz[k] += a*(fz[k]+ Fz_t) ; 
	/*   	printf("from extforce_vel_Up k=%d name=%s vz =%f Fz_t=%f\n", k, gp->name,vz[k],Fz_t); */
}
void extforce_parms(GROUP *gp)
{
	EXTFORCE_PARMS *parm;
	char *string; 
	parm = ddcCalloc(1, sizeof(EXTFORCE_PARMS));
	gp->itype = EXTFORCE;
        object_get((OBJECT *) gp, "Fzeq", &string, STRING, 1, NULL);
        parm->Fzeq = eq_parse(string,"m*l/t^2","t");
	gp->parm = parm; 
	gp->write_dynamics = NULL;
	gp->velocityUpdate= (void (*)(int,int,GROUP*,void*,double,double))extforce_velocityUpdate; 
	gp->Update= (void (*) (GROUP*, int, void *, double, double))extforce_Update; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
