#include "nglfNK.h"
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
#include "preduce.h"
#include "group.h"
#include "gid.h"
#include "codata.h"
#include "units.h"
#include "three_algebra.h"


void kinetic_terms(SYSTEM*sys, int flag);
void eval_energyInfo(SYSTEM *sys);

NGLFNK_PARMS *nglfNK_parms(INTEGRATOR*integrator)
{
   NGLFNK_PARMS *parms;
   parms = ddcMalloc(sizeof(NGLFNK_PARMS));
   object_get((OBJECT*)integrator,"P",  &parms->P,  WITH_UNITS,1,"0.0","pressure",NULL); 
   object_get((OBJECT*)integrator,"W",  &parms->W,  WITH_UNITS,3,"1.0 1.0 1.0","mass",NULL); 
   object_get((OBJECT*)integrator,"T",  &parms->T,  WITH_UNITS,1,"1.0","T",NULL); 
   object_get((OBJECT*)integrator,"tau",&parms->tau,WITH_UNITS,1,"1.0","t",NULL); 
   return parms;
}
void nglfNK(DDC*ddc, SIMULATE*simulate, NGLFNK_PARMS *p)
{
   double dt = simulate->dt;
   double halfdt = 0.5*dt; 
   double time = simulate->time;
   SYSTEM* sys = simulate->system;
   STATE* state = sys->collection->state;
   unsigned nlocal= state->nlocal; 
   BOX_STRUCT* box =sys->box; 

   THREE_SMATRIX *sion = &(sys->energyInfo.sion);
   //THREE_SMATRIX *tion = &(sys->energyInfo.tion);
   THREE_VECTOR L,dLdt;
   THREE_VECTOR S[nlocal],dSdt[nlocal];
   THREE_MATRIX h=box_get_h(sys->box); 
   THREE_MATRIX dhdt=box_get_dhdt(sys->box); 
   THREE_VECTOR corner = box_get_corner(box); 
   THREE_VECTOR reducedCorner = box_get_reducedcorner(box); 

/*    Get canonical variables   */
   VSET(L,h.xx,h.yy,h.zz);
   VSET(dLdt,dhdt.xx,dhdt.yy,dhdt.zz);

   double V = L.x*L.y*L.z; 
   for (unsigned kk=0; kk<nlocal; kk++) 
   {
      S[kk].x = fmod((state->rx[kk]-corner.x)/L.x,1.0); 
      S[kk].y = fmod((state->ry[kk]-corner.y)/L.y,1.0); 
      S[kk].z = fmod((state->rz[kk]-corner.z)/L.z,1.0); 
      dSdt[kk].x = (state->vx[kk] - state->rx[kk]*dLdt.x/L.x)/L.x;
      dSdt[kk].y = (state->vy[kk] - state->ry[kk]*dLdt.y/L.y)/L.y;
      dSdt[kk].z = (state->vz[kk] - state->rz[kk]*dLdt.z/L.z)/L.z;
   }

/* Integrate canonical velocities to have half timestep */
   THREE_VECTOR P; 
   VSET(P,-(sion->xx),-(sion->yy),-(sion->zz)); 
// double cP = units_convert(1.0,NULL,"bar"); 
   P.x = P.y = 0.5*(P.x+P.y); 
   THREE_SMATRIX tion0 = (sys->energyInfo.tion);
   double kBT = kB*p->T; 
   double mu =1.0/p->tau; 
   //mu = 0.0; 
   RANDOM *random=system_getRandom(sys);
   for (unsigned kk=0; kk<nlocal; kk++) 
   {
      void *randomParms = random_getParms(random, kk); 
      THREE_VECTOR g = gasdev3d(random,randomParms); 
      double rmass = 1.0/((ATOMTYPE_PARMS *) (state->species[kk]->parm))->mass;
      double sigma = sqrt(2.0*kBT*(rmass*mu)/halfdt);
      dSdt[kk].x = dSdt[kk].x + halfdt*((state->fx[kk]*rmass-mu*dLdt.x*S[kk].x + sigma*g.x)-( mu*L.x + 2*dLdt.x)*dSdt[kk].x)/L.x;
      dSdt[kk].y = dSdt[kk].y + halfdt*((state->fy[kk]*rmass-mu*dLdt.y*S[kk].y + sigma*g.y)-( mu*L.y + 2*dLdt.y)*dSdt[kk].y)/L.y;
      dSdt[kk].z = dSdt[kk].z + halfdt*((state->fz[kk]*rmass-mu*dLdt.z*S[kk].z + sigma*g.z)-( mu*L.z + 2*dLdt.z)*dSdt[kk].z)/L.z;
   }
   dLdt.x += halfdt*V/(p->W.x*L.x)*(P.x - p->P); 
   dLdt.y += halfdt*V/(p->W.y*L.y)*(P.y - p->P); 
   dLdt.z += halfdt*V/(p->W.z*L.z)*(P.z - p->P); 

/* integrate canical positions to full timestep */
   for (unsigned kk = 0; kk < nlocal; kk++)
   {
      S[kk].x =  fmod(S[kk].x + dt*dSdt[kk].x,1.0);
      S[kk].y =  fmod(S[kk].y + dt*dSdt[kk].y,1.0);
      S[kk].z =  fmod(S[kk].z + dt*dSdt[kk].z,1.0);
   }
   h.xx=L.x +=  dt*dLdt.x;
   h.yy=L.y +=  dt*dLdt.y;
   h.zz=L.z +=  dt*dLdt.z;
   V = L.x*L.y*L.z; 

/*  Transform from canonical positions  to native positions */
   box_put(box,HO,&h); 
   for (unsigned kk =0;kk<nlocal;kk++) 
   {
      state->rx[kk] = L.x*(S[kk].x+reducedCorner.x); 
      state->ry[kk] = L.y*(S[kk].y+reducedCorner.y); 
      state->rz[kk] = L.z*(S[kk].z+reducedCorner.z); 
   }
   time += dt;                                     // positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt 
   simulate->time=sys->time=time; 
   simulate->loop++;
   sys->loop = simulate->loop;
   ddc->update = 0;

   if (ddcenergy(ddc, sys, 0) != 0) return;

   /* Integrate canonical velocities to full  timestep */
   nlocal=state->nlocal; 
   VSET(P,-(sion->xx-tion0.xx)/V,-(sion->yy-tion0.yy)/V,-(sion->zz-tion0.zz)/V);           // Volume
   P.x = P.y = 0.5*(P.x+P.y); 

   dLdt.x +=    halfdt*V/(p->W.x*L.x)*(P.x - p->P); 
   dLdt.y +=    halfdt*V/(p->W.y*L.y)*(P.y - p->P); 
   dLdt.z +=    halfdt*V/(p->W.z*L.z)*(P.z - p->P); 

   for (unsigned kk=0; kk<nlocal; kk++) 
   {
      void *randomParms = random_getParms(random, kk); 
      THREE_VECTOR g = gasdev3d(random,randomParms); 
      double rmass = 1.0/((ATOMTYPE_PARMS *) (state->species[kk]->parm))->mass;
      //double sigma = sqrt(2.0*(dt/2)*kBT*(rmass*mu));
      double sigma = sqrt(2.0*kBT*(rmass*mu)/halfdt);
      dSdt[kk].x = (dSdt[kk].x + halfdt*(state->fx[kk]*rmass-mu*dLdt.x*S[kk].x + sigma*g.x)/L.x)/(1.0+halfdt*( mu*L.x + 2*dLdt.x)/L.x);
      dSdt[kk].y = (dSdt[kk].y + halfdt*(state->fy[kk]*rmass-mu*dLdt.y*S[kk].y + sigma*g.y)/L.y)/(1.0+halfdt*( mu*L.y + 2*dLdt.y)/L.y);
      dSdt[kk].z = (dSdt[kk].z + halfdt*(state->fz[kk]*rmass-mu*dLdt.z*S[kk].z + sigma*g.z)/L.z)/(1.0+halfdt*( mu*L.z + 2*dLdt.z)/L.z);

   }

   /*  Transform from canonical velocities  to native velocities */
   dhdt = mzero; 
   dhdt.xx=dLdt.x;
   dhdt.yy=dLdt.y;
   dhdt.zz=dLdt.z;
   box_put(box,DHDT,&dhdt); 
   for (unsigned kk =0;kk<nlocal;kk++) 
   {
      state->vx[kk] = L.x*dSdt[kk].x+S[kk].x*dLdt.x; 
      state->vy[kk] = L.y*dSdt[kk].y+S[kk].y*dLdt.y; 
      state->vz[kk] = L.z*dSdt[kk].z+S[kk].z*dLdt.z; 
   }
   kinetic_terms(sys, 1);
   eval_energyInfo(sys);

}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
