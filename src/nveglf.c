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
#include "simulate.h"
#include "ddcMalloc.h"
#include "ddcenergy.h"
#include "preduce.h"

typedef struct nveglf_parms_st
{
	void *ptr;
} NVEGLF_PARMS;

void ddckinetic(int mode, SYSTEM*sys);

NVEGLF_PARMS *nveglf_parms(INTEGRATOR*integrator)
{
	NVEGLF_PARMS *parms;
	parms = ddcMalloc(sizeof(NVEGLF_PARMS));
	return parms;
}

void nveglf(DDC*ddc, SIMULATE*simulate, NVEGLF_PARMS*p)
{
	double time, dt, fac, rmass, vol_atom, volume;
	unsigned nlocal;
	gid_type nsimul;
	STATE *state;
	SYSTEM *sys;
	double *rx, *ry, *rz;
	double *vx, *vy, *vz;
	double *fx, *fy, *fz;
	sys = simulate->system;
	nsimul = sys->nglobal;
	nlocal = sys->nlocal;
	vol_atom = sys->box->volume/nsimul;
	state = sys->collection->state;
	dt = simulate->dt;
	time = simulate->time;
	rx = state->rx;
	ry = state->ry;
	rz = state->rz;
	vx = state->vx;
	vy = state->vy;
	vz = state->vz;
	fx = state->fx;
	fy = state->fy;
	fz = state->fz;
	rmass = 1.0/((ATOMTYPE_PARMS *) (sys->species[0]->parm))->mass;
	if (simulate->Veq != NULL)
	{
		vol_atom = (simulate->Veq->function) (time, (void *)simulate->Veq);
		volume = vol_atom*sys->nglobal;
		fac = cbrt(volume/sys->box->volume);
		for (unsigned kk = 0; kk < nlocal; kk++)
		{
			rx[kk] *= fac;
			ry[kk] *= fac;
			rz[kk] *= fac;
		}
		box_put(sys->box, VOLUME, &volume);
	}
		for (unsigned kk = 0; kk < nlocal; kk++)
		{		/* [L_3a] */
			vx[kk] += 0.5*rmass*dt*fx[kk];
			vy[kk] += 0.5*rmass*dt*fy[kk];
			vz[kk] += 0.5*rmass*dt*fz[kk];
		}
		for (unsigned kk = 0; kk < nlocal; kk++)
		{		/* [L_0] */
			rx[kk] += dt*vx[kk];
			ry[kk] += dt*vy[kk];
			rz[kk] += dt*vz[kk];
		}
		if (simulate->Veq != NULL)
		{
			vol_atom = (simulate->Veq->function) (time + dt, (void *)simulate->Veq);
			volume = vol_atom*sys->nglobal;
			fac = cbrt(volume/sys->box->volume);
			for (unsigned kk = 0; kk < nlocal; kk++)	/* [L_5][L_4] */
			{
				rx[kk] *= fac;
				ry[kk] *= fac;
				rz[kk] *= fac;
			}
			box_put(sys->box, VOLUME, &volume);
		}
		for (unsigned kk = 0; kk < nlocal; kk++)
		   backInBox(rx + kk, ry + kk, rz + kk);
		ddc->update = 0;
		time += dt;
		simulate->loop++;
		sys->loop = simulate->loop;
		ddcenergy(ddc, sys, 1);
/*     positions, vol_atom, box (volume, h0,hinv ) , and forces at  t = n*dt + dt */
		nlocal = sys->nlocal;
		rx = state->rx;
		ry = state->ry;
		rz = state->rz;
		vx = state->vx;
		vy = state->vy;
		vz = state->vz;
		fx = state->fx;
		fy = state->fy;
		fz = state->fz;

		for (unsigned kk = 0; kk < nlocal; kk++)
		{		/* [L_3a] */
			vx[kk] += 0.5*rmass*simulate->dt*fx[kk];
			vy[kk] += 0.5*rmass*simulate->dt*fy[kk];
			vz[kk] += 0.5*rmass*simulate->dt*fz[kk];
		}
		ddckinetic(0, sys);
		/*errorCheck(ddc->domain_id, simulate->loop, state, e, p, datafile); */
	simulate->time = sys->time = time;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
