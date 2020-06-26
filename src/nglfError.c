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
#include "units.h"

void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);

NGLF_PARMS *nglfError_parms(INTEGRATOR*integrator)
{
	NGLF_PARMS *parms;
	parms = ddcMalloc(sizeof(NGLF_PARMS));

	return parms;
}
void timestep(DDC*ddc, SIMULATE*simulate, NGLF_PARMS*p, double dt)
{
	SYSTEM *sys=simulate->system;
	STATE *state=sys->collection->state; 
   double time = sys->time; 
	ddcenergy(ddc, sys, -1);
	for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update(sys->group[kk],FRONT_TIMESTEP,state,sys->time,0.5*dt);
	for (unsigned kk=0; kk<sys->nlocal; kk++) 
	{
		GROUP* group = state->group[kk]; 
		group->velocityUpdate(FRONT_TIMESTEP,kk,group,state,time,0.5*dt);
	}
	for (unsigned kk = 0; kk < sys->nlocal; kk++)
	{
	   state->rx[kk] += dt*state->vx[kk];
	   state->ry[kk] += dt*state->vy[kk];
	   state->rz[kk] += dt*state->vz[kk];
	}
	for (unsigned kk = 0; kk < sys->nlocal; kk++) backInBox_fast(state->rx + kk, state->ry + kk, state->rz + kk);
	ddc->update = 0;
	sys->time+=dt; 
	ddcenergy(ddc, sys, -1);
	for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update(sys->group[kk],BACK_TIMESTEP,state,time,0.5*dt);
	for (unsigned kk=0;kk<sys->nlocal;kk++) 
	{
		GROUP* group = state->group[kk]; 
		group->velocityUpdate(BACK_TIMESTEP,kk,group,state,time,0.5*dt);
	}
	kinetic_terms(sys, 1);
	eval_energyInfo(sys);
}
void saveState0(SYSTEM *sys,int *nlocal,double *time, THREE_VECTOR *r,THREE_VECTOR *v, THREE_VECTOR *f,double *K, double *U)
{
	STATE* state = sys->collection->state;
	*time = sys->time; 
	*nlocal = sys->nlocal; 
	for (unsigned i=0;i<sys->nlocal;i++)
	{
		VSET(r[i],state->rx[i],state->ry[i],state->rz[i]); 
		VSET(v[i],state->vx[i],state->vy[i],state->vz[i]); 
		VSET(f[i],state->fx[i],state->fy[i],state->fz[i]); 
		double mass = ((ATOMTYPE_PARMS *) (state->species[i]->parm))->mass;
		double Ki = 0.5*mass*(SQ(state->vx[i])+SQ(state->vy[i])+SQ(state->vz[i])); 
		double Ui = state->potentialEnergy[i];  
		K[i] = Ki; 
		U[i] = Ui; 
	}
}
void restoreState0(SYSTEM *sys,int *nlocal,double *time, THREE_VECTOR *r,THREE_VECTOR *v, THREE_VECTOR *f,double *K, double *U)
{
	STATE* state = sys->collection->state;
	for (unsigned i=0;i<sys->nlocal;i++)
	{
		state->rx[i]=r[i].x;state->ry[i]=r[i].y;state->rz[i]=r[i].z; 
		state->vx[i]=v[i].x;state->vy[i]=v[i].y;state->vz[i]=v[i].z; 
		state->fx[i]=f[i].x;state->fy[i]=f[i].y;state->fz[i]=f[i].z; 
	}
	sys->time = *time; 
	sys->nlocal=*nlocal; 
}
void errorTimeStep(SYSTEM *sys,double *K0, double *U0)
{
		int nlocal = sys->nlocal; 
		STATE*state = sys->collection->state; 
		double errU=0.0; 
		double errK=0.0; 
		double errUmax=0.0; 
		double errKmax=0.0; 
	double energy_convert = units_convert(1.0,NULL,"eV");
		double energy=0.0; 
		double energy0=0.0; 
		double errTotalMax=0.0,errTotalMin=0.0; 
		double errTotalBar=0.0; 
		double errTotalBarSmall=0.0; 
		for (int i=0;i<nlocal;i++) 
		{
			double mass = ((ATOMTYPE_PARMS *) (state->species[i]->parm))->mass;
			double Ki = 0.5*mass*(SQ(state->vx[i])+SQ(state->vy[i])+SQ(state->vz[i])); 
			double Ui = state->potentialEnergy[i];  
			double dKi = (K0[i]-Ki)*energy_convert; 
			double dUi = (U0[i]-Ui)*energy_convert;
			errK += dKi  ; 
			errU += dUi  ; 
			if (errKmax < fabs(dKi) )  errKmax = fabs(dKi); 
			if (errUmax < fabs(dUi) )  errUmax = fabs(dUi); 
			double errTotal = dKi + dUi; 
			if (errTotalMin > errTotal || i==0) errTotalMin = errTotal;  
			if (errTotalMax < errTotal || i==0) errTotalMax = errTotal;  
			errTotalBar += errTotal; 
			energy0 += U0[i]  ; 
			energy+=Ui; 
			//printf("%d %e %e %e %e\n",i,dKi,dUi,errKmax,errUmax); 
		}
		errTotalBar /= nlocal; 
		errU /= nlocal; 
		errK /= nlocal; 
		static double delta,emin ;
		static FILE *file=NULL,*file1=NULL,*file2=NULL; 
		static double bin[2][100]; 
		static double escale=0.0; 
		if (file == NULL) 
		{
			escale  = 0.7*MAX(fabs(errTotalMax),fabs(errTotalMin)); 
			delta = 1.5*(errTotalMax-errTotalMin)/100.0; 
			emin =  1.5*errTotalMin; 
			for (int i=0;i<100;i++)  bin[0][i]=bin[1][i]=0.0 ;
			file=fopen("error.data","w"); 
		}
		if (file1== NULL) file1=fopen("time.data","w"); 
		if (file2== NULL) file2=fopen("data.data","w"); 
		static int cnt=0; 
		for (int i=0;i<nlocal;i++) 
		{
			double mass = ((ATOMTYPE_PARMS *) (state->species[i]->parm))->mass;
			double ai = sqrt(SQ(state->fx[i])+SQ(state->fy[i])+SQ(state->fz[i]))/mass; 
			double Ki = 0.5*mass*(SQ(state->vx[i])+SQ(state->vy[i])+SQ(state->vz[i])); 
			double Ui = state->potentialEnergy[i];  
			double dKi = (K0[i]-Ki)*energy_convert; 
			double dUi = (U0[i]-Ui)*energy_convert;
			double errTotal = dKi + dUi; 
		    int ibin = (errTotal - emin)/delta; 
			if (ibin < 0 ) ibin =0; 
			if (ibin > 99) ibin =99; 
			bin[0][ibin]+=1.0 ; 
			if ( fabs(errTotal)< escale) errTotalBarSmall += errTotal; 
			
			fprintf(file2,"%e %e %e %e %e %e\n",errTotal,dKi,dUi,Ki,Ui,ai); 
		}
		errTotalBarSmall /= nlocal; 
		errU /= nlocal; 
		cnt++; 
  		fseek(file, 0, SEEK_SET);
		for (int i=0;i<100;i++)  
		{
			double e = emin + delta *(i+0.5); 
			double p = bin[0][i]/(100*delta*cnt) ;  
			fprintf(file,"%f %e %e\n",e,p,p*e); 
		}
	fprintf(file,"end_of_data\n"); 
	fprintf(file2,"end_of_data\n"); 
	fprintf(file1,"%"PRId64" %f %e %e %e %e",sys->loop,sys->time,errTotalBar,errTotalBarSmall,errK,errU); fprintf(file1," %e %e\n",errKmax,errUmax); 
	fflush(file); 
	fflush(file1); 
	fflush(file2); 
}


void nglfError(DDC*ddc, SIMULATE*simulate, NGLF_PARMS*p)
{
	SYSTEM* sys = simulate->system;
	double dt = simulate->dt;
	SIGNED64 loop = simulate->loop;

   	THREE_VECTOR r0[sys->nlocal],v0[sys->nlocal],f0[sys->nlocal]; 
	double t0,K0[sys->nlocal],U0[sys->nlocal]; 
	int nlocal0; 

   	THREE_VECTOR r1[sys->nlocal],v1[sys->nlocal],f1[sys->nlocal]; 
	double t1,K1[sys->nlocal],U1[sys->nlocal]; 
	int nlocal1; 

	ddcenergy(ddc, sys, 0);
	saveState0(sys,&nlocal0,&t0,r0,v0,f0,K0,U0); 
	loop = simulate->loop; 

	timestep(ddc, simulate, p,dt);
	saveState0(sys,&nlocal1,&t1,r1,v1,f1,K1,U1); 

	restoreState0(sys,&nlocal0,&t0,r0,v0,f0,K0,U0); 
	sys->loop = loop; 
	simulate->loop=loop;  
	simulate->time=t0;  
	sys->time=t0;  

	for (int l=0;l<2;l++) timestep(ddc, simulate, p,0.5*dt);

	sys->loop++; 
	simulate->time=sys->time; 
	simulate->loop=sys->loop; 

	errorTimeStep(sys,K1,U1); 
}

