#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "three_algebra.h"
#include "particle.h"
#include "system.h"
#include "units.h"
#include "neighbor.h"
#include "ptiming.h"

void mgpt(SYSTEM*sys, NBR*nbr, void *, ETYPE*e);
void kinetic_terms(SYSTEM*sys, int flag);
void mgptparmVol(double vol);

/*
void energy(SYSTEM*sys, GEOM*geom)
{
	static THREE_VECTOR vzero = { 0.0, 0.0, 0.0 };
	double vol, voln;
	static PARTICLESET *particleset = NULL;
	static NBR *nbr;
	STATE *state;
	ETYPE *e = &sys->e; 
    unsigned i, nion, nlocal;
    gid_type nsimul;
	nion = sys->nion;
	nlocal = sys->nlocal;
	nsimul = sys->nglobal;
	state = sys->collection->state;
	vol = sys->box->volume;
	voln = vol/nion;
	mgptparmVol(voln);

        particleset = ParticleSet(particleset, state->rx, state->ry, state->rz, 8, state->label, 8, NULL, 4, &nion, &nion);
	for (i = 0; i < nion; i++)
	   backInBox(state->rx + i, state->ry + i, state->rz + i);
	nbr = neighbors(geom, particleset);
	for (i = 0; i < nbr->npairs; i++) nbr->pairs[i].p->fp = vzero;
	mgpt(sys, nbr, sys->potential[0]->parms, e);
	kinetic_terms(sys, 1);
	e->pion += (2.0*e->rk)/(3.0);
	e->pion /= vol;
	SMATNORM(e->sion, vol);
}
*/

/* Add kinetic contributions to stress tensor and energy If flag <1, then only compute kinetic energy */
void kinetic_terms(SYSTEM*sys, int flag)
{
   profile(KINETIC_TERMS,START);
	ETYPE *e = &sys->energyInfo; 
	ETYPE *ge;
	ETYPE *se;
	THREE_VECTOR J = { 0.0,0.0,0.0}; 
	THREE_SMATRIX *sion = sys->collection->state->sion;
	double *potential = sys->collection->state->potentialEnergy; 
	double *vx = sys->collection->state->vx;
	double *vy = sys->collection->state->vy;
	double *vz = sys->collection->state->vz;
	unsigned nlocal = sys->nlocal;
	gid_type nsimul = sys->nglobal;

	e->rk=0.0; 
	e->number=0.0; 
   e->tion.xx =0.0;       
   e->tion.yy =0.0;
   e->tion.zz =0.0;
   e->tion.xy =0.0;
   e->tion.xz =0.0;
   e->tion.yz =0.0;
   for (int i=0;i<sys->ngroup;i++)
   {
      ge = &sys->group[i]->energyInfo;
      ge->tion.xx =0.0;       
      ge->tion.yy =0.0;
      ge->tion.zz =0.0;
      ge->tion.xy =0.0;
      ge->tion.xz =0.0;
      ge->tion.yz =0.0;
   }
   for (int i=0;i<sys->nspecies;i++)
   {
      se = &sys->species[i]->energyInfo;
      se->tion.xx =0.0;       
      se->tion.yy =0.0;
      se->tion.zz =0.0;
      se->tion.xy =0.0;
      se->tion.xz =0.0;
      se->tion.yz =0.0;
   }

   for (unsigned k = 0; k < nlocal; k++)
   {
      double mass = ((ATOMTYPE_PARMS *) (sys->collection->state->species[k]->parm))->mass;
      double vxx = vx[k]*vx[k]; 
      double vyy = vy[k]*vy[k]; 
      double vzz = vz[k]*vz[k]; 
      double vxy = vx[k]*vy[k]; 
      double vxz = vx[k]*vz[k]; 
      double vyz = vy[k]*vz[k]; 

      double K = 0.5*mass*(vxx + vyy + vzz);
      double U = potential[k]; 
      THREE_SMATRIX S = sion[k];
      J.x += (K + U)*vx[k]-0.5*(S.xx*vx[k] + S.xy*vy[k] + S.xz*vz[k]) ;
      J.y += (K + U)*vy[k]-0.5*(S.xy*vx[k] + S.yy*vy[k] + S.yz*vz[k]) ;
      J.z += (K + U)*vz[k]-0.5*(S.xz*vx[k] + S.yz*vy[k] + S.zz*vz[k]) ;
      e->tion.xx += mass*vxx;
      e->tion.yy += mass*vyy;
      e->tion.zz += mass*vzz;
      e->tion.xy += mass*vxy;
      e->tion.xz += mass*vxz;
      e->tion.yz += mass*vyz;
      sion[k].xx -= mass*vxx;
      sion[k].yy -= mass*vyy;
      sion[k].zz -= mass*vzz;
      sion[k].xy -= mass*vxy;
      sion[k].xz -= mass*vxz;
      sion[k].yz -= mass*vyz;
      e->mass += mass; 
      e->number++; 
      e->rk += K;


      ge = &sys->collection->state->group[k]->energyInfo;
      ge->tion.xx += mass*vxx;
      ge->tion.yy += mass*vyy;
      ge->tion.zz += mass*vzz;
      ge->tion.xy += mass*vxy;
      ge->tion.xz += mass*vxz;
      ge->tion.yz += mass*vyz;
      ge->mass += mass; 
      ge->number++; 
      ge->rk += K;
      ge->eion += U;

      se = &sys->collection->state->species[k]->energyInfo;
      se->tion.xx += mass*vxx;
      se->tion.yy += mass*vyy;
      se->tion.zz += mass*vzz;
      se->tion.xy += mass*vxy;
      se->tion.xz += mass*vxz;
      se->tion.yz += mass*vyz;
      se->mass += mass; 
      se->number++; 
      se->rk += K;
      se->eion += U;
   }
   e->thermal_flux = J; 

   e->temperature = 2.0*e->rk/(3.0*nsimul); 
   for (int i=0;i<sys->ngroup;i++)
   {
      ge = &sys->group[i]->energyInfo;
      /*ge->temperature = 2.0*ge->rk/(3.0*nsimul);*/
   }
   for (int i=0;i<sys->nspecies;i++)
   {
      se = &sys->species[i]->energyInfo;
      /*ge->temperature = 2.0*ge->rk/(3.0*nsimul);*/
   }
   profile(KINETIC_TERMS,END);
}

void kinetic_scaleCorrection(ETYPE*e, double scale)
{
   double s2;
   s2 = scale*scale;
   e->rk *= s2;
   e->temperature *= s2;
   e->sion.xx = e->sion.xx + (s2 - 1.0)*e->tion.xx;
   e->sion.yy = e->sion.yy + (s2 - 1.0)*e->tion.yy;
   e->sion.zz = e->sion.zz + (s2 - 1.0)*e->tion.zz;
   e->sion.xy = e->sion.xy + (s2 - 1.0)*e->tion.xy;
   e->sion.xz = e->sion.xz + (s2 - 1.0)*e->tion.xz;
   e->sion.yz = e->sion.yz + (s2 - 1.0)*e->tion.yz;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
