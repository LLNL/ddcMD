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

void	VPreduceOrthorhombicB7_OneLatticeReduction(double *rx, double *ry, double *rz, unsigned nlocal);
void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);
NGLF_PARMS *nveglf_simple_parms(INTEGRATOR*integrator)
{
	NGLF_PARMS *parms;
	parms = ddcMalloc(sizeof(NGLF_PARMS));
	return parms;
}
void nveglf_simple(DDC*ddc, SIMULATE*simulate, NGLF_PARMS*p)
{
	double time, dt;
	unsigned nlocal;
	STATE *state;
	SYSTEM *sys;
	double mass; 
	SPECIES **species; 
	double *rx, *ry, *rz;
	double *vx, *vy, *vz;
	double *fx, *fy, *fz;
	sys = simulate->system;
	nlocal = sys->nlocal;
	state = sys->collection->state;
	dt = simulate->dt;
	time = simulate->time;
	rx = state->rx; ry = state->ry; rz = state->rz;
	vx = state->vx; vy = state->vy; vz = state->vz;
	fx = state->fx; fy = state->fy; fz = state->fz;
	species = state->species; 
	double hdt = 0.5*dt; 
/*      All Variables at t = n*dt;  */
		for (unsigned kk=0;kk<sys->nlocal;kk++) 
		{
     		mass = ((ATOMTYPE_PARMS *) (species[kk]->parm))->mass;
			double a = hdt/mass; 
		   vx[kk] += a*fx[kk];
		   vy[kk] += a*fy[kk];
		   vz[kk] += a*fz[kk];
		}
		for (unsigned kk = 0; kk < nlocal; kk++)
		{
		   rx[kk] += dt*vx[kk];
		   ry[kk] += dt*vy[kk];
		   rz[kk] += dt*vz[kk];
		}
		ddc->update = 0;
/*     positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt */
		time += dt;
		simulate->time=sys->time=time; 
		simulate->loop++;
		sys->loop = simulate->loop;
		if (ddcenergy(ddc, sys, 0) != 0) return;

		nlocal = sys->nlocal;
		rx = state->rx; ry = state->ry; rz = state->rz;
		vx = state->vx; vy = state->vy; vz = state->vz;
		fx = state->fx; fy = state->fy; fz = state->fz;
		species = state->species; 

		for (unsigned kk=0;kk<sys->nlocal;kk++) 
		{
			mass = ((ATOMTYPE_PARMS *) (species[kk]->parm))->mass;
			double a = hdt/mass; 
		   vx[kk] += a*fx[kk];
		   vy[kk] += a*fy[kk];
		   vz[kk] += a*fz[kk];
		}
	simulate->time = sys->time = time;
	//for (unsigned kk = 0; kk < nlocal; kk++) Preduce_fast(rx + kk, ry + kk, rz + kk);
	kinetic_terms(sys, 1);
	eval_energyInfo(sys);
	VPreduceOrthorhombicB7_OneLatticeReduction(rx, ry, rz, nlocal);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
