#include "nptglf.h"
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
#include "langevin.h"
#include "ddcMalloc.h"
#include "ddcenergy.h"
#include "lcg64.h"
#include "units.h"
#include "collection.h"
#include "preduce.h"

int getRank(int);
void ddckinetic(int mode, SYSTEM*sys);
void kinetic_scaleCorrection(ETYPE *, double);
NPTGLF_PARMS *nptglf_parms(INTEGRATOR*integrator)
{
	NPTGLF_PARMS *parms;
	parms = ddcMalloc(sizeof(NPTGLF_PARMS));
	object_get((OBJECT *) integrator, "Gamma", &parms->Gamma, WITH_UNITS, 1, "1","m/l^4",NULL);
	object_get((OBJECT *) integrator, "zeta", &parms->zeta, WITH_UNITS, 1, "1","pressure*t",NULL);
	object_get((OBJECT *) integrator, "pressure", &parms->pressure, WITH_UNITS, 1, "1","pressure",NULL);
	return parms;
}

void nptglf_writedynamic(INTEGRATOR *integrator, FILE *file)
{
		double zeta = ((NPTGLF_PARMS *) integrator->parms)->zeta;
		fprintf(file, "%s INTEGRATOR { zeta=%16.12e ; }\n", integrator->name, units_convert(zeta,NULL,"pressure*t"));
}

void nptglf(DDC*ddc, SIMULATE*simulate, NPTGLF_PARMS*p)
{
	double time, dt, fac, vol_atom, volume;
	double Gamma, zeta, Peq, deltap;
	double *rx, *ry, *rz;
	double *vx, *vy, *vz;
	SYSTEM *sys = simulate->system;
	ETYPE *energyInfo = &sys->energyInfo;
	STATE *state = sys->collection->state;
	gid_type nsimul = sys->nglobal;
	unsigned nlocal = sys->nlocal;
	vol_atom = sys->box->volume/nsimul;
	dt = simulate->dt;
	time = simulate->time;
	Gamma = p->Gamma;
	zeta = p->zeta;
	rx = state->rx; ry = state->ry; rz = state->rz;
	vx = state->vx; vy = state->vy; vz = state->vz;
	//Peq = (simulate->Peq->function) (time, (void *)simulate->Peq);
	Peq = p->pressure;
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
		deltap = energyInfo->pion - Peq;
		zeta = (zeta + 0.5*dt*deltap);
		fac = exp(-zeta*dt/(6.0*Gamma*vol_atom));
		for (unsigned kk = 0; kk < nlocal; kk++)
		{
			vx[kk] = vx[kk]*fac ;
			vy[kk] = vy[kk]*fac ;
			vz[kk] = vz[kk]*fac ;
		}
		for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update(sys->group[kk],FRONT_TIMESTEP,state,time,0.5*dt);
		for (unsigned kk = 0; kk < nlocal; kk++)
		{
			GROUP *group = state->group[kk]; 
			group->velocityUpdate(FRONT_TIMESTEP,kk,group, state, time,0.5*dt);
		}
		vol_atom += 0.5*dt/Gamma*zeta;
		fac = exp(zeta*dt/(6.0*Gamma*vol_atom));
//      velocities , zeta and vol_atom at t = n*dt + 0.5d ;  
		for (unsigned kk = 0; kk < nlocal; kk++)
		{		
			rx[kk] = (fac*rx[kk] +dt*vx[kk])*fac;
			ry[kk] = (fac*ry[kk] +dt*vy[kk])*fac;
			rz[kk] = (fac*rz[kk] +dt*vz[kk])*fac;
		}
		vol_atom += 0.5*dt/Gamma*zeta;
		volume = vol_atom*nsimul;
		box_put(sys->box, VOLUME, &volume); 
		for (unsigned kk = 0; kk < nlocal; kk++) backInBox(rx + kk, ry + kk, rz + kk);
		ddc->update = 0;
		time += dt;
		simulate->time = sys->time = time;
		simulate->loop++;
		sys->loop = simulate->loop;
		ddcenergy(ddc, sys, 1);
		nlocal = sys->nlocal; 
		//Peq = (simulate->Peq->function) (time, (void *)simulate->Peq);
      Peq=p->pressure;
		deltap = energyInfo->pion - Peq;
//     positions, vol_atom, box volume, h0, & hinv  , and forces at  t = n*dt + dt 
		rx = state->rx; ry = state->ry; rz = state->rz;
		vx = state->vx; vy = state->vy; vz = state->vz;
		for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update(sys->group[kk],BACK_TIMESTEP,state,time,0.5*dt);
		for (unsigned kk = 0; kk < nlocal; kk++)
		{
			GROUP *group = state->group[kk]; 
			group->velocityUpdate(BACK_TIMESTEP,kk,group, state, time,0.5*dt);
		}
		{
			double rk, zeta0, p0, fac_old;
			rk = energyInfo->rk;
			ddckinetic(0, sys);
			energyInfo->pion = energyInfo->pion + 2.0/3.0*(energyInfo->rk - rk)/volume;
			p0 = energyInfo->pion;
			deltap = energyInfo->pion - Peq;
			zeta0 = zeta;
			zeta = (zeta0 + 0.5*dt*deltap);
			fac = exp(- zeta*dt/(6.0*Gamma*vol_atom));
			rk = energyInfo->rk;
			for (unsigned kk = 0; kk < 5; kk++)
			{
				energyInfo->pion = p0 + (fac*fac - 1.0)*(2.0/3.0)*rk/volume;
				deltap = energyInfo->pion - Peq;
				zeta = (zeta0 + 0.5*dt*deltap);
				fac_old = fac;
				fac = exp(-zeta*dt/(6.0*Gamma*vol_atom));
				if (fabs(fac_old/fac - 1.0) < 1e-12) break;
			}
		}
		for (unsigned kk = 0; kk < nlocal; kk++)	// [L_5][L_4] 
		{
			vx[kk] *= fac;
			vy[kk] *= fac;
			vz[kk] *= fac;
		}
		kinetic_scaleCorrection(energyInfo, fac);

		// FILE *datafile = stdout; 
		//errorCheck(ddc->domain_id, simulate->loop, state, e, p, datafile); 
	simulate->time = sys->time = time;
	p->zeta = zeta;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
