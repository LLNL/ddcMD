#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "system.h"
#include "ddcMalloc.h"
#include "units.h"
#include "preduce.h"
#include "external.h"


typedef struct slice_st
{
	double center, width, set_velocity, set_temp;
	double temperature,chi,delta,chi_backward,delta_backward;
	double  rk,vf,af,num_atoms, mass;
	THREE_VECTOR v,f,v_backward;
} SLICE;

typedef struct shear_parms_st
{
	double tau;
	SLICE *top, *bottom;
} SHEAR_PARMS;
char *shear_printinfo(GROUP *g,FILE *file)
{
	static char line[1024];
	double  dTtop,dTbot,dVtop,dVbot; 
	extern int getRank(int);
	SHEAR_PARMS *p; 
	SYSTEM *sys  ;
	sys = g->parent; 
	if (sys->loop % 100 ) return NULL; 
	
	SLICE *ptop,*pbot;
	p=(SHEAR_PARMS *)g->parm ;
	ptop = p->top;
	pbot = p->bottom;
	dTtop = units_convert(ptop->temperature-ptop->set_temp,NULL,"T"); 
	dVtop = units_convert(ptop->v.y-ptop->set_velocity,NULL,"velocity"); 
	dTbot = units_convert(pbot->temperature-pbot->set_temp,NULL,"T"); 
	dVbot = units_convert(pbot->v.y-pbot->set_velocity,NULL,"velocity"); 
	sprintf(line,"%"PRId64" %f %f %f %f",sys->loop,dTtop,dVtop,dTbot,dVbot); 
	if (getRank(0) == 0) {fputs(line,file); printf("\n"); }
	return line; 
}


void shearreduce(SHEAR_PARMS *e)
{
	double psum[22], sum[22];
	int i; 
	i=0; 
	psum[i] = e->top->v.x; i++; 
	psum[i] = e->top->v.y; i++; 
	psum[i] = e->top->v.z; i++; 
	psum[i] = e->top->f.x; i++; 
	psum[i] = e->top->f.y; i++; 
	psum[i] = e->top->f.z; i++; 
	psum[i] = e->top->vf; i++; 
	psum[i] = e->top->af; i++; 
	psum[i] = e->top->rk; i++; 
	psum[i] = e->top->num_atoms; i++; 
	psum[i] = e->top->mass; i++; 

	psum[i] = e->bottom->v.x; i++; 
	psum[i] = e->bottom->v.y; i++; 
	psum[i] = e->bottom->v.z; i++; 
	psum[i] = e->bottom->f.x; i++; 
	psum[i] = e->bottom->f.y; i++; 
	psum[i] = e->bottom->f.z; i++; 
	psum[i] = e->bottom->vf; i++; 
	psum[i] = e->bottom->af; i++; 
	psum[i] = e->bottom->rk; i++; 
	psum[i] = e->bottom->num_atoms; i++; 
	psum[i] = e->bottom->mass; i++; 
	MPI_Allreduce(psum, sum, i, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
	i=0; 
	e->top->v.x = sum[i++];
	e->top->v.y = sum[i++];
	e->top->v.z = sum[i++];
	e->top->f.x = sum[i++];
	e->top->f.y = sum[i++];
	e->top->f.z = sum[i++];
	e->top->vf = sum[i++];
	e->top->af = sum[i++];
	e->top->rk = sum[i++];
	e->top->num_atoms = sum[i++];
	e->top->mass = sum[i++];

	e->bottom->v.x = sum[i++];
	e->bottom->v.y = sum[i++];
	e->bottom->v.z = sum[i++];
	e->bottom->f.x = sum[i++];
	e->bottom->f.y = sum[i++];
	e->bottom->f.z = sum[i++];
	e->bottom->vf = sum[i++];
	e->bottom->af = sum[i++];
	e->bottom->rk = sum[i++];
	e->bottom->num_atoms = sum[i++];
	e->bottom->mass = sum[i++];
}
void shear_Update(GROUP *g, int mode, STATE *state, double time, double dt)
{
   double mass,chi,c,temperature,x,y,ztop,zbot;
   int k, nlocal;
   SHEAR_PARMS *p; 
   SLICE *ptop, *pbot;
   p=(SHEAR_PARMS *)g->parm ;
   ptop = p->top;
   pbot = p->bottom;
   double *rz = state->rz; 
   double *vx = state->vx; double *vy = state->vy; double *vz = state->vz; 
   double *fx = state->fx; double *fy = state->fy; double *fz = state->fz; 
   rz = state->rz;
   nlocal = state->nlocal;
   ptop->af=ptop->vf=ptop->mass=ptop->num_atoms=ptop->rk=0.0;
   ptop->v = vzero; 
   ptop->f = vzero; 


   pbot->af=pbot->vf=pbot->mass=pbot->num_atoms=pbot->rk=0.0;
   pbot->v = vzero; 
   pbot->f = vzero; 
   x=y=0.0; 

   for (k=0;k<nlocal;k++)
   {
      ztop = rz[k]-p->top->center; backInBox(&x,&y,&ztop); 
      zbot = rz[k]-p->bottom->center; backInBox(&x,&y,&zbot); 
      mass = ((ATOMTYPE_PARMS *)(state->species[k]->parm))->mass; 
      if(fabs(ztop) < 0.5*ptop->width)
      {
         ptop->num_atoms += 1.0;
         ptop->mass += mass; 
         ptop->v.x += mass*vx[k];
         ptop->v.y += mass*vy[k];
         ptop->v.z += mass*vz[k];
         ptop->f.x += fx[k];
         ptop->f.y += fy[k];
         ptop->f.z += fz[k];
         ptop->af += (fx[k]*fx[k]+fy[k]*fy[k]+fz[k]*fz[k])/mass;
         ptop->vf += vx[k]*fx[k] + vy[k]*fy[k] + vz[k]*fz[k];
         ptop->rk +=0.5*mass*(vx[k]*vx[k]+vy[k]*vy[k]+vz[k]*vz[k]);

      }
      else if(fabs(zbot) < 0.5*pbot->width)
      {
         pbot->num_atoms+=1.0;
         pbot->mass += mass; 
         pbot->v.x += mass*vx[k];
         pbot->v.y += mass*vy[k];
         pbot->v.z += mass*vz[k];
         pbot->f.x += fx[k];
         pbot->f.y += fy[k];
         pbot->f.z += fz[k];
         pbot->af += (fx[k]*fx[k]+fy[k]*fy[k]+fz[k]*fz[k])/mass;
         pbot->vf += vx[k]*fx[k] + vy[k]*fy[k] + vz[k]*fz[k];
         pbot->rk +=0.5*mass*(vx[k]*vx[k]+vy[k]*vy[k]+vz[k]*vz[k]);
      }
   }
   shearreduce(p);

   ptop->v.x /=ptop->mass;
   ptop->v.y /=ptop->mass;
   ptop->v.z /=ptop->mass;
   ptop->rk -= 0.5*ptop->mass*VSQ(ptop->v);
   ptop->temperature = 2.0*ptop->rk/(3.0*(ptop->num_atoms-1.0)); 
   ptop->delta = dt/p->tau * (ptop->set_velocity - ptop->v.y);
   ptop->chi = sqrt(1.0+ (dt/p->tau)*(ptop->set_temp/ptop->temperature -1.0));

   ptop->v_backward.x = ( ptop->v.x + dt*ptop->f.x /ptop->mass );
   ptop->v_backward.y = ( ptop->v.y + dt*ptop->f.y /ptop->mass + (dt/p->tau) * ptop->set_velocity )/(1.0+dt/p->tau); 
   ptop->v_backward.z = ( ptop->v.z + dt*ptop->f.z /ptop->mass );
   ptop->delta_backward = dt/p->tau * (ptop->set_velocity - ptop->v_backward.y);
   c = ptop->temperature + ( 2*dt*( ptop->vf - DOT(ptop->v,ptop->f) ) + dt*dt*( ptop ->af - VSQ(ptop->f)/ptop->mass))/(3.0*ptop->num_atoms-3.0); 

   temperature = ptop->temperature; 
   for (k=0;k<5;k++) 
   {
      chi = sqrt(1.0 + (dt/p->tau)*(ptop->set_temp/temperature -1.0));
      temperature  = c/((2-chi)*(2.0-chi));
      /*printf("temperature=%f\n",temperature);  */
   }
   ptop->chi_backward = chi ; 

   pbot->v.x /=pbot->mass;
   pbot->v.y /=pbot->mass;
   pbot->v.z /=pbot->mass;
   pbot->rk -= 0.5*pbot->mass*VSQ(pbot->v);
   pbot->temperature = 2.0*pbot->rk/(3.0*(pbot->num_atoms-1)); 

   pbot->delta = dt/p->tau * (pbot->set_velocity - pbot->v.y);
   pbot->chi = sqrt(1.0+ (dt/p->tau)*(pbot->set_temp/pbot->temperature -1.0));

   pbot->v_backward.x = ( pbot->v.x + dt*pbot->f.x /pbot->mass );
   pbot->v_backward.y = ( pbot->v.y + dt*pbot->f.y /pbot->mass + (dt/p->tau) * pbot->set_velocity )/(1.0+dt/p->tau); 
   pbot->v_backward.z = ( pbot->v.z + dt*pbot->f.z /pbot->mass );
   pbot->delta_backward = dt/p->tau * (pbot->set_velocity - pbot->v_backward.y);
   c = pbot->temperature + ( 2*dt*( pbot->vf - DOT(pbot->v,pbot->f) ) + dt*dt*( pbot ->af - VSQ(pbot->f)/pbot->mass))/(3.0*pbot->num_atoms-3.0); 

   temperature = pbot->temperature; 
   for (k=0;k<5;k++) 
   {
      chi = sqrt(1.0 + (dt/p->tau)*(pbot->set_temp/temperature -1.0));
      temperature  = c/((2-chi)*(2.0-chi));
      /*printf("temperature=%f\n",temperature);  */
   }
   pbot->chi_backward = chi ; 
   shear_printinfo(g,stdout); 
}
void shear_velocityUpdate(int mode, int k,GROUP *g, STATE *state, double time, double dt)
{
   double *rz = state->rz; 
   double *vx = state->vx; 
   double *vy = state->vy; 
   double *vz = state->vz; 
   double *fx = state->fx; 
   double *fy = state->fy; 
   double *fz = state->fz; 
   SPECIES **species = state->species; 
   double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;
   double chi,delta,rmass;
   double x,y,ztop,zbot; 
   THREE_VECTOR v; 
   SHEAR_PARMS *p; 
   p=g->parm ;
   rmass = 1.0/mass; 
   chi=1.0; 
   delta =0; 
   v = vzero; 
   x=y=0.0; 
   ztop = rz[k]-p->top->center; backInBox(&x,&y,&ztop); 
   zbot = rz[k]-p->bottom->center; backInBox(&x,&y,&zbot); 
   if (mode==FRONT_TIMESTEP)
   {
      if(fabs(ztop) < 0.5*p->top->width)
      {
         v = p->top->v; 
         chi = p->top->chi; 
         delta = p->top->delta; 

      }
      if(fabs(zbot) < 0.5*p->bottom->width)
      {
         v = p->bottom->v; 
         chi = p->bottom->chi; 
         delta = p->bottom->delta; 
      }
      vx[k] +=  dt*fx[k]*rmass;
      vy[k] +=  dt*fy[k]*rmass;
      vz[k] +=  dt*fz[k]*rmass;
      vx[k] +=  (chi-1.0)*(vx[k]-v.x);
      vy[k] +=  (chi-1.0)*(vy[k]-v.y)+delta;
      vz[k] +=  (chi-1.0)*(vz[k]-v.z);
   }
   if (mode==BACK_TIMESTEP)
   {
      if(fabs(ztop) < 0.5*p->top->width)
      {
         v = p->top->v_backward; 
         chi = p->top->chi_backward; 
         delta = p->top->delta_backward; 
      }
      if(fabs(zbot) < 0.5*p->bottom->width)
      {
         v = p->bottom->v_backward; 
         chi = p->bottom->chi_backward; 
         delta = p->bottom->delta_backward; 
      }
      vx[k] +=  dt*fx[k]*rmass;
      vy[k] +=  dt*fy[k]*rmass;
      vz[k] +=  dt*fz[k]*rmass;
      vx[k] +=  (chi-1.0)*(vx[k]-v.x);
      vy[k] +=  (chi-1.0)*(vy[k]-v.y)+delta;
      vz[k] +=  (chi-1.0)*(vz[k]-v.z);
   }
}
void shear_parms(GROUP *gp)
{
   SHEAR_PARMS *parm;
   SLICE *ptop, *pbot;
   parm = ddcCalloc(1, sizeof(SHEAR_PARMS));
   ptop = ddcCalloc(1, sizeof(SLICE));
   pbot = ddcCalloc(1, sizeof(SLICE));
   gp->itype = SHEAR;
   object_get((OBJECT *) gp, "tau", &parm->tau, WITH_UNITS, 1, "1.0","t",NULL);
   object_get((OBJECT *) gp, "top_center", &ptop->center , WITH_UNITS, 1, "-1","l",NULL);
   object_get((OBJECT *) gp, "bottom_center", &pbot->center , WITH_UNITS, 1, "-1","l",NULL);
   object_get((OBJECT *) gp, "top_width", &ptop->width , WITH_UNITS, 1, "-1","l",NULL);
   object_get((OBJECT *) gp, "bottom_width", &pbot->width , WITH_UNITS, 1, "-1","l",NULL);
   object_get((OBJECT *) gp, "top_velocity", &ptop->set_velocity , WITH_UNITS, 1, "-1","l/t",NULL);
   object_get((OBJECT *) gp, "bottom_velocity", &pbot->set_velocity , WITH_UNITS, 1, "-1","l/t",NULL);
   object_get((OBJECT *) gp, "top_temp", &ptop->set_temp , WITH_UNITS, 1, "-1","T",NULL);
   object_get((OBJECT *) gp, "bottom_temp", &pbot->set_temp , WITH_UNITS, 1, "-1","T",NULL);


   parm->top=ptop;
   parm->bottom=pbot;
   gp->parm = parm;
   gp->write_dynamics = NULL;
   gp->velocityUpdate= (void (*)(int, int, GROUP*, void *, double, double))shear_velocityUpdate; 
   gp->Update= (void (*) (GROUP*, int, void *, double, double))shear_Update; 
   /*
      printf("shear parms: %f %f %f %f\n",ptop->center,ptop->width,pbot->center,pbot->width);
      printf("more parms : %f %f %f %f\n",ptop->set_velocity,ptop->set_temp,pbot->set_velocity,pbot->set_temp);
      */

}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
