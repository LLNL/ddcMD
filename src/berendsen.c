#include "berendsen.h"
#include <math.h>
#include <string.h>
#include "state.h"
#include "ddcMalloc.h"
#include "object.h"
#include "system.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))

typedef struct Berendsen_parms_st
{
   double Teq;
   double tau;
   unsigned interval;
   double lambda;
   double Tsum;  // sum of temperature for all steps over an interval
   unsigned nT;  // number of temperatures in the sum
   unsigned doScaling;
} BERENDSEN_PARMS;



/** Calculates the velocity scaling factor parms->lambda to use in the
 *  next call of velocityUpdate (which will be at the FRONT_TIMESTEP of
 *  the next step).  Note that Tave is either the average temperature
 *  over some interval (when parms->interval > 1) or equivalent to the
 *  instantaneous temperature (when parms->interval==1).
 */
void berendsen_Update(GROUP* g, int mode, void* state, double time, double dt)
{
   
   if (mode == FRONT_TIMESTEP)
   {
      BERENDSEN_PARMS* parms = g->parm;
      parms->Tsum += g->energyInfo.temperature;
      parms->nT += 1;
      double Tave = parms->Tsum/parms->nT;
      double tau = parms->tau;
      double Teq = parms->Teq;

      double ratio; // ratio = Teq/Tave, but Teq or Tave might be zero.
      if (Tave == 0)
	 ratio = 0; // Arbitrary value.  Scaling v=0 => v=0.
      else
	 ratio = Teq/Tave;

      // note that dt here is 1/2 of actual time step.
      if (tau != 0)
	 parms->lambda = sqrt(1 + (2.0*dt/tau) * ((ratio)-1));
      else
	 parms->lambda = sqrt(ratio);
      
      parms->doScaling = 0;
      if (system_getLoop(NULL) % parms->interval == 0)
      {
	 parms->Tsum = 0;
	 parms->nT = 0;
	 parms->doScaling = 1;
      }
   }
}

void berendsen_velocityUpdate(int mode, int k, GROUP* g, STATE *state, double time, double dt)
{
   double *vx = state->vx; 
   double *vy = state->vy; 
   double *vz = state->vz; 
   double *fx = state->fx; 
   double *fy = state->fy; 
   double *fz = state->fz; 
   SPECIES **species = state->species; 
   double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;
   if (mode == FRONT_TIMESTEP)
   {
      BERENDSEN_PARMS* parms = g->parm;
      if (parms->doScaling == 1)
      {
	 vx[k] *= parms->lambda;
	 vy[k] *= parms->lambda;
	 vz[k] *= parms->lambda;
      }
   }
   
   double a = dt/mass; 
   vx[k] += a*fx[k];
   vy[k] += a*fy[k];
   vz[k] += a*fz[k];
}

void berendsen_parms(GROUP* gp)
{
   BERENDSEN_PARMS* parms = ddcMalloc(sizeof(BERENDSEN_PARMS));

   OBJECT* obj = (OBJECT*) gp;

   object_get(obj, "Teq", &parms->Teq, WITH_UNITS, 1, "0.0", "T", NULL);
   char* tmp;
   object_get(obj, "tau", &tmp, STRING, 1, "foo");
   if (strcmp(tmp, "dt") == 0)
      parms->tau = 0.0;
   else
      object_get(obj, "tau", &parms->tau, WITH_UNITS, 1, "1.0", "t", NULL);
   object_get(obj, "interval", &parms->interval, INT, 1, "1");
   ddcFree(tmp);
   
   gp->itype = BERENDSEN;
   gp->parm = parms;
   gp->velocityUpdate = (void (*)(int, int, GROUP *, void *, double, double))berendsen_velocityUpdate; 
   gp->Update = berendsen_Update; 
   parms->nT = 0;
   parms->Tsum = 0.0;
   parms->doScaling = 0;
}
