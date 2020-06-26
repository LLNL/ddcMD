#include "nglf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "ddc.h"
#include "box.h"
#include "species.h"
#include "ddcMalloc.h"
#include "ddcenergy.h"
#include "expandbuffer.h"
#include "auxNeighbor.h"
#include "preduce.h"
#include "group.h"

void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);
void nglf_collision(SYSTEM *sys, double dt);
double pair_function(SYSTEM *sys, int i,int j,double r, double *dr);
typedef struct  ref_str {THREE_VECTOR r,v;} REF ; 
REF *vaf_v0(void); 
static AUXNEIGHBOR  *auxNeighbor ; 
NGLF_PARMS *nglf_parms(INTEGRATOR*integrator)
{
	NGLF_PARMS *parms;
	parms = ddcMalloc(sizeof(NGLF_PARMS));
    double rmax = 3.0;   
    auxNeighbor = auxNeighbor_request(rmax);  
	return parms;
}
void scalePositionsByBoxChange(BOX_STRUCT *box, double time, double *rx, double *ry, double *rz, unsigned nlocal)
{
	THREE_MATRIX hfac; 
	THREE_VECTOR rold,r; 
	box_put(box,BOX_TIME,(void *)&time);
	box_get(box,HFAC,(void *)&hfac);
	if ( ! matrix_equal(hfac, I_3x3))
	{
     		for (unsigned kk = 0; kk < nlocal; kk++)
     		{
			rold.x=rx[kk];
			rold.y=ry[kk];
			rold.z=rz[kk];
			r = matrix_vector(hfac,rold); 
			rx[kk]=r.x;
			ry[kk]=r.y;
			rz[kk]=r.z;
		
     		}
	}
}
void	updateStateAliases(SYSTEM *sys,unsigned *nlocal, double **rx,double **ry,double **rz,double **vx,double **vy,double **vz,double **fx,double **fy,double **fz,SPECIES ***species,gid_type **label)
{
	STATE *state=sys->collection->state; 
	*nlocal = sys->nlocal;
	*rx = state->rx; *ry = state->ry; *rz = state->rz; // The SYSTEM and STATE might change during the call to ddcenergy
	*vx = state->vx; *vy = state->vy; *vz = state->vz; // (i.e., we might reassign particles to new tasks) so we need to
	*fx = state->fx; *fy = state->fy; *fz = state->fz; // update all of the aliases we use.
	*label = state->label; 
	*species = state->species; 
}
void nglf(DDC*ddc, SIMULATE*simulate, NGLF_PARMS*p)
{
	double dt = simulate->dt;
	double time = simulate->time;
	SYSTEM* sys = simulate->system;
	STATE* state = sys->collection->state;
   scalePositionsByBoxChange(sys->box,time,state->rx,state->ry,state->rz,state->nlocal); 
   for (unsigned kk=0; kk<sys->nlocal; kk++) 
   {
      GROUP* group = state->group[kk]; 
      group->velocityUpdate(FRONT_TIMESTEP,kk,group,state,time,0.5*dt);
   }
   REF *v0 = vaf_v0(); 
   for (int  kk = 0; kk < state->nlocal; kk++)
   {
      THREE_VECTOR delta; 
      state->rx[kk] += delta.x = dt*state->vx[kk];
      state->ry[kk] += delta.y = dt*state->vy[kk];
      state->rz[kk] += delta.z = dt*state->vz[kk];
      if (v0 != NULL )  VOP1(v0[kk].r,+=,delta); 
   }
   double time_plus_dt = time + dt; 
   scalePositionsByBoxChange(sys->box,time_plus_dt,state->rx,state->ry,state->rz,state->nlocal); 
   for (int kk = 0; kk < state->nlocal; kk++) backInBox_fast(state->rx + kk, state->ry + kk, state->rz + kk);
   ddc->update = 0;
   time += dt;                                     // positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt 
   simulate->time=sys->time=time; 
   simulate->loop++;
   sys->loop = simulate->loop;
   for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update1(sys->group[kk],-1,state,time,0.5*dt);
   if (ddcenergy(ddc, sys, 0) != 0) return;
   for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update(sys->group[kk],BACK_TIMESTEP,state,time,0.5*dt);

   for (unsigned kk=0;kk<sys->nlocal;kk++) 
   {
      GROUP* group = state->group[kk]; 
      group->velocityUpdate(BACK_TIMESTEP,kk,group,state,time,0.50*dt);
   }
   kinetic_terms(sys, 1);
   for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update2(sys->group[kk],-1,state,time,0.5*dt);
   //eval_energyInfo(sys);
   for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update(sys->group[kk],FRONT_TIMESTEP,state,time,0.5*dt);

   /*errorCheck(ddc->domain_id, simulate->loop, state, sys->energyInfo, p, datafile); */
   simulate->time = sys->time = time;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
