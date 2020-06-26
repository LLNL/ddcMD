#include "nglfTest.h"
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
#include "codata.h"
#include "preduce.h"

#define NBIN 50 

typedef struct istate_str {THREE_VECTOR r,v,f,dfdt;double K,U;} ISTATE; 
void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);
void step(DDC *ddc, SYSTEM *sys, double dt);
void stepm(DDC *ddc, SYSTEM *sys, double dt,int nInner);
void saveIstate(ISTATE *state,unsigned n);
void loadIstate(ISTATE *state,unsigned n);
void errorIstate(ISTATE *state0, ISTATE *state_h,ISTATE *state,ISTATE *ref,unsigned n, char *distName,int nerr,double bin[NBIN][2]);

static double *rx, *ry, *rz;
static double *vx, *vy, *vz;
static double *fx, *fy, *fz;
static double *q; 
static double *energy; 
static SPECIES **species; 
static unsigned nlocal ;
static double dt; 
static AUXNEIGHBOR  *auxNeighbor ; 
static AUXNEIGHBORLIST *nlist = NULL;
static AUXNEIGHBORINDEX *nindex = NULL;
static NGLFTEST_PARMS *_parms ;
static double binM[NBIN][2]; 
static  int nerrM =0; 
static double binS[NBIN][2]; 
static  int nerrS =0; 
NGLFTEST_PARMS *nglfTest_parms(INTEGRATOR*integrator)
{
	NGLFTEST_PARMS *parms=NULL;
	parms = ddcMalloc(sizeof(NGLFTEST_PARMS));
	parms->acc = NULL; 
   object_get((OBJECT *) integrator, "subDivide", &parms->subDivide, INT, 1, "1");
   object_get((OBJECT *) integrator, "highAccuarcyDt", &parms->highAccuarcyDt, WITH_UNITS, 1, "0.25","t",NULL);
   object_get((OBJECT *) integrator, "singleDt", &parms->singleDt, WITH_UNITS, 1, "1.00","t",NULL);
   object_get((OBJECT *) integrator, "rmax", &parms->rmax, WITH_UNITS, 1, "0.0","l",NULL);
	double rmax = 2.0;
    auxNeighbor = auxNeighbor_request(rmax); 
	_parms = parms; 

	return parms;
}

void nglfTest(DDC*ddc, SIMULATE*simulate, NGLFTEST_PARMS*parms)
{
	SYSTEM *sys = simulate->system;
	STATE *state = sys->collection->state;
	species = state->species; 
	nlocal = sys->nlocal;
	unsigned nion = sys->nion;
	dt = simulate->dt;
	int nI=  dt /parms->highAccuarcyDt + .9999999; 
	int nSingle=  dt /parms->singleDt +  .9999999; 
	double Dt = dt/nI; 
	double time = simulate->time;
	ISTATE state0[nion],state1[nion],state_Mstep[nion],state_Sstep[nion],state_ref[nion]; 
	rx = state->rx; ry = state->ry; rz = state->rz;
	vx = state->vx; vy = state->vy; vz = state->vz;
	fx = state->fx; fy = state->fy; fz = state->fz;
	q= state->q; 
	energy = state->potentialEnergy; 
	_parms = parms; 
	for (int loop = 0; loop < 1; loop++) // was loop over printrate
	{
		saveIstate(state0,nion); 
		stepm(ddc,sys,dt,parms->subDivide); 
		saveIstate(state_Mstep,nion); 

		loadIstate(state0,nion); 
		for (int kk=0;kk<nSingle ;kk++) step(ddc,sys,dt/nSingle); 
		saveIstate(state_Sstep,nion); 

		loadIstate(state0,nion); 
	    for (int kk=0;kk<nI;kk++) step(ddc,sys,0.5*Dt); 
		saveIstate(state1,nion); 
	    for (int kk=0;kk<nI;kk++) step(ddc,sys,0.5*Dt); 
		saveIstate(state_ref,nion); 

		errorIstate(state0,state1,state_Sstep,state_ref,nlocal,"SingleStep.dist",nerrS,binS); 
		errorIstate(state0,state1,state_Mstep,state_ref,nlocal,"MultiStep.dist",nerrM,binM); 
		time += dt;
		simulate->time=sys->time=time; 
		simulate->loop++;
		sys->loop = simulate->loop;
	}
	for (unsigned  kk = 0; kk < nlocal; kk++) backInBox_fast(rx + kk, ry + kk, rz + kk);
	ddc->update = 0;
	ddcenergy(ddc, sys, 0);
	kinetic_terms(sys, 1);
	eval_energyInfo(sys);
}
void step(DDC *ddc, SYSTEM *sys, double dt)
{
	for (unsigned kk = 0; kk < nlocal; kk++)
	{
	   double	a = 0.5*dt/((ATOMTYPE_PARMS *) (species[kk]->parm))->mass;
	   rx[kk]+=dt*(vx[kk]+a*fx[kk]);
	   ry[kk]+=dt*(vy[kk]+a*fy[kk]);
	   rz[kk]+=dt*(vz[kk]+a*fz[kk]);
	   vx[kk]+=a*fx[kk];
	   vy[kk]+=a*fy[kk];
	   vz[kk]+=a*fz[kk];
	}
	ddcenergy(ddc, sys, -1);
	for (unsigned kk = 0; kk < nlocal; kk++)
	{
	   double	a = 0.5*dt/((ATOMTYPE_PARMS *) (species[kk]->parm))->mass;
	   vx[kk]+=a*fx[kk];
	   vy[kk]+=a*fy[kk];
	   vz[kk]+=a*fz[kk];
	}
}
THREE_VECTOR forceFast(int i)
{
	THREE_VECTOR ff = { 0.0,0.0,0.0}; 
    int start = nindex[i].startPairs; 
    int np = nindex[i].nPairs; 
	double rmax = _parms->rmax;
	double rmaxI = 1.0/rmax; 
	
	for (int kk=0;kk<np;kk++)
	{
		int j = (nlist[start+kk].j); 
		double dx = rx[i] - rx[j];
		double dy = ry[i] - ry[j];
		double dz = rz[i] - rz[j];
		double r = sqrt(dx*dx+dy*dy+dz*dz); 
		nearestImage(&dx,&dy,&dz); 
		double sr = r*rmaxI; 
		double s = (sr < 1.0) ? (1.0-sr*sr*sr)*(1.0-sr*sr*sr): 0.0; 
		ff.x += ke*q[i]*q[j]/(r*r*r) * s* dx ; 
		ff.y += ke*q[i]*q[j]/(r*r*r) * s* dy ; 
		ff.z += ke*q[i]*q[j]/(r*r*r) * s* dz ; 
	}
	return ff; 
}
void stepm(DDC *ddc, SYSTEM *sys, double dt,int nInner)
{
	THREE_VECTOR ff[nlocal]; 
	nlist = auxNeighbor_list();
	nindex =auxNeighbor_index();
	for (unsigned i = 0; i < nlocal; i++) ff[i]= forceFast(i);
	for (unsigned i = 0; i < nlocal; i++)
	{
	   double	a = 0.5*dt/((ATOMTYPE_PARMS *) (species[i]->parm))->mass;
	   vx[i]+=a*(fx[i]-ff[i].x);
	   vy[i]+=a*(fy[i]-ff[i].y);
	   vz[i]+=a*(fz[i]-ff[i].z);
	}
	for (int kk=0;kk<nInner;kk++)
	{
		for (unsigned i = 0; i < nlocal; i++)
		{
	   		double	a = 0.5*dt/(nInner*((ATOMTYPE_PARMS *) (species[i]->parm))->mass);
	   		vx[i]+=a*(ff[i].x);
	   		vy[i]+=a*(ff[i].y);
	   		vz[i]+=a*(ff[i].z);
		}
		for (unsigned i = 0; i < nlocal; i++)
		{
			double a=dt/nInner; 
	   		rx[i]+=a*(vx[i]);
	   		ry[i]+=a*(vy[i]);
	   		rz[i]+=a*(vz[i]);
		}
		for (unsigned i = 0; i < nlocal; i++) ff[i]= forceFast(i);
		for (unsigned i = 0; i < nlocal; i++)
		{
	   		double	a = 0.5*dt/(nInner*((ATOMTYPE_PARMS *) (species[i]->parm))->mass);
	   		vx[i]+=a*(ff[i].x);
	   		vy[i]+=a*(ff[i].y);
	   		vz[i]+=a*(ff[i].z);
		}
	}
	ddcenergy(ddc, sys, -1);
	for (unsigned i = 0; i < nlocal; i++)
	{
	   double	a = 0.5*dt/((ATOMTYPE_PARMS *) (species[i]->parm))->mass;
	   vx[i]+=a*(fx[i]-ff[i].x);
	   vy[i]+=a*(fy[i]-ff[i].y);
	   vz[i]+=a*(fz[i]-ff[i].z);
	}
}
void errorIstate(ISTATE *state0, ISTATE *state_h, ISTATE *state,ISTATE *ref,unsigned n,char *distName,int nerr, double bin[NBIN][2])
{
	static FILE *file=NULL; 
	if (file == NULL) file = fopen("error.dat","w"); 
	double dUmax=0.0; 
	double dKmax=0.0; 
	double dUbar=0.0; 
	double dKbar=0.0; 
	THREE_VECTOR df,ddf;
	AUXNEIGHBORLIST *nlist = auxNeighbor_list();
	AUXNEIGHBORINDEX *nindex =auxNeighbor_index();
	THREE_VECTOR fs; 
	if (nerr==0) for (int i=0;i<50;i++) bin[i][0]=bin[i][1]=0.0; 
	for (unsigned i=0;i<n;i++) 
	{
		int start = nindex[i].startPairs; 
		int j = (nlist[start].j); 
		double r = sqrt(nlist[start].r2); 
		double dx = ref[i].r.x - ref[j].r.x;
		double dy = ref[i].r.y - ref[j].r.y;
		double dz = ref[i].r.z - ref[j].r.z;
		double r1= r; 
		nearestImage(&dx,&dy,&dz); 
		fs.x = ref[i].f.x - ke*q[i]*q[j]/(r*r*r) * dx ; 
		fs.y = ref[i].f.y - ke*q[i]*q[j]/(r*r*r) * dy ; 
		fs.z = ref[i].f.z - ke*q[i]*q[j]/(r*r*r) * dz ; 
	    double	m = ((ATOMTYPE_PARMS *) (species[i]->parm))->mass;
		double dU= fabs(state[i].U-ref[i].U); 
		double dK= fabs(state[i].K-ref[i].K); 
		//double f = sqrt(VSQ(ref[i].f)); 
		double v = sqrt(VSQ(ref[i].v)); 
		df.x= ref[i].f.x-state0[i].f.x; 
		df.y= ref[i].f.y-state0[i].f.y; 
		df.z= ref[i].f.z-state0[i].f.z; 
		ddf.x= ref[i].f.x+state0[i].f.x-2.0*state_h[i].f.x; 
		ddf.y= ref[i].f.y+state0[i].f.y-2.0*state_h[i].f.y; 
		ddf.z= ref[i].f.z+state0[i].f.z-2.0*state_h[i].f.z; 
		double Df = sqrt(VSQ(df))/(dt*m);
		double DDf = sqrt(VSQ(ddf))/(dt*dt*m);
		dUmax = MAX(dUmax,dU); 
		dKmax = MAX(dKmax,dK); 
		dUbar += dU; 
		dKbar += dK; 
		fprintf(file,"%d %f %f %f %f %f %f %e %e\n",i,v,r1,1/(r1*r1),sqrt(VSQ(fs)),Df,DDf,fabs(dK),fabs(dU));
		fflush(file); 
		int ierr = -5.0*log(dK)/log(10.); 
		ierr = MIN(50,ierr); 
		ierr = MAX(0,ierr); 
		bin[ierr][0]+=1.0; 
		bin[ierr][1]+=dK; 
		nerr++; 
	}
	dKbar /= n; 
	dUbar /= n; 
	printf("error max =%e %e %e %e\n",dKmax,dUmax,dKbar,dUbar); 
	FILE * distfile = fopen(distName,"w"); 
	for (int i=0;i<50;i++) fprintf(distfile,"%f %lf %e\n",-i*0.2,bin[i][0]/nerr,bin[i][1]/nerr); 
	fclose(distfile); 
//fclose(file); 
}
void saveIstate(ISTATE *state,unsigned n)
{	
	for (unsigned i=0;i<n;i++) 
	{
	   double	mass = ((ATOMTYPE_PARMS *) (species[i]->parm))->mass;
		state[i].r.x=rx[i]; 
		state[i].r.y=ry[i]; 
		state[i].r.z=rz[i]; 
		state[i].v.x=vx[i]; 
		state[i].v.y=vy[i]; 
		state[i].v.z=vz[i]; 
		state[i].f.x=fx[i]; 
		state[i].f.y=fy[i]; 
		state[i].f.z=fz[i]; 
		state[i].U=energy[i]; 
		state[i].K=0.5*mass*VSQ(state[i].v); 
	}
}
void loadIstate(ISTATE *state,unsigned n)
{	
	for (unsigned i=0;i<n;i++) 
	{
		rx[i]=state[i].r.x; 
		ry[i]=state[i].r.y; 
		rz[i]=state[i].r.z; 
		vx[i]=state[i].v.x; 
		vy[i]=state[i].v.y; 
		vz[i]=state[i].v.z; 
		fx[i]=state[i].f.x; 
		fy[i]=state[i].f.y; 
		fz[i]=state[i].f.z; 
		energy[i]=state[i].U; 
	}
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
