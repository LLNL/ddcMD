#include "addVelocity.h"

#include <assert.h>
#include <mpi.h>

#include "object.h"
#include "species.h"
#include "group.h"
#include "simulate.h"
#include "ddcMalloc.h"
#include "units.h"
#include "three_algebra.h"

int getRank(int);
static THREE_VECTOR computeVcm(SYSTEM* sys);


typedef struct addVelocity_parms_st
{
   int by_species;
   double* vx;
   double* vy;
   double* vz;
   
} ADD_VELOCITY_PARMS;

void* addVelocity_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   ADD_VELOCITY_PARMS* parms = ddcMalloc(sizeof(ADD_VELOCITY_PARMS));
   parms->by_species = -1;
   
   char** speciesNames = NULL;
   char** groupNames = NULL;
   int nsIn = object_getv(obj, "species", (void*)&speciesNames, STRING,IGNORE_IF_NOT_FOUND);
   int ngIn = object_getv(obj, "groups",  (void*)&groupNames,   STRING,IGNORE_IF_NOT_FOUND);

   if (nsIn > 0 && ngIn > 0 )
   {
      if (getRank(0) == 0)
	 printf("Error in addvelocity transform object.\n"
		"Please specify either groups or species, not both.\n");
      exit(2);
   }

   int nv;
   int velSize = 0;
   if (nsIn > 0)
   {
      nv=system_getNspecies(NULL);
      velSize = 3*nsIn;
      parms->by_species = 1;
   }
   if (ngIn > 0)
   {
      group_get(NULL, NGROUPS, (void*) &nv);
      velSize = 3*ngIn;
      parms->by_species = 0;
   }
   
   if (parms->by_species < 0 )
   {
      if (getRank(0) == 0)
	 printf("Error in addvelocity transform object.\n"
		"Please specify either groups or species.\n");
      exit(2);
   }
   
   parms->vx = ddcCalloc(nv, sizeof(double));
   parms->vy = ddcCalloc(nv, sizeof(double));
   parms->vz = ddcCalloc(nv, sizeof(double));
   for (int ii=0; ii<nv; ++ii)
   {
      parms->vx[ii] = 0;
      parms->vy[ii] = 0;
      parms->vz[ii] = 0;
   }
   
   double* vel = NULL;
   int nVel = object_getv(obj, "velocity", (void*)&vel, DOUBLE,IGNORE_IF_NOT_FOUND);
   if (nVel != velSize)
   {
      if (getRank(0) == 0)
	 printf("Error in addvelocity transform object.\n"
		"Please specify 3 velocity components for "
		"each group or species\n");
      exit(2);
   }
   
   double length_convert = units_convert(1.0,"l",NULL); 
   double time_convert = units_convert(1.0,"t",NULL); 
   double velocity_convert = length_convert/time_convert; 
   for (int ii=0; ii<nsIn; ++ii)
   {
      SPECIES* species = species_find(NULL, speciesNames[ii]);
      if (species)
      {
	 int iSpecies = species->index;
	 parms->vx[iSpecies] = vel[3*ii] * velocity_convert;
	 parms->vy[iSpecies] = vel[3*ii + 1] * velocity_convert;
	 parms->vz[iSpecies] = vel[3*ii + 2] * velocity_convert;
      }
      ddcFree(speciesNames[ii]);
   }
   for (int ii=0; ii<ngIn; ++ii)
   {
      GROUP* group = group_find(NULL, groupNames[ii]);
      if (group)
      {
	 int iGroup = group->index;
	 parms->vx[iGroup] = vel[3*ii] * velocity_convert;
	 parms->vy[iGroup] = vel[3*ii + 1] * velocity_convert;
	 parms->vz[iGroup] = vel[3*ii + 2] * velocity_convert;
      }
      ddcFree(groupNames[ii]);
   }

   ddcFree(speciesNames);
   ddcFree(groupNames);
   ddcFree(vel);
   return parms;
}

/** You can only set vcm for the system as a whole, not for individual
 * groups or species. */
void* setVelocity_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system;
   
   ADD_VELOCITY_PARMS* parms = ddcMalloc(sizeof(ADD_VELOCITY_PARMS));
   parms->by_species = -1; 

   double* vTarget = (double*) ddcMalloc(3*sizeof(double));
   int nvIn = object_get(obj, "vcm", (void*)vTarget, WITH_UNITS, 3, "0.0 0.0 0.0", "velocity", NULL);

   if (nvIn != 3 )
   {
      if (getRank(0) == 0)
	 printf("Error in setvelocity transform object.\n"
		"Please specify 3 components for the center of mass velocity\n");
      exit(2);
   }

   int nv = system_getNspecies(sys);
   parms->by_species = 1;
   
   parms->vx = ddcCalloc(nv, sizeof(double));
   parms->vy = ddcCalloc(nv, sizeof(double));
   parms->vz = ddcCalloc(nv, sizeof(double));
   for (int ii=0; ii<nv; ++ii)
   {
      parms->vx[ii] = 0;
      parms->vy[ii] = 0;
      parms->vz[ii] = 0;
   }

   THREE_VECTOR vcm = computeVcm(sys);
   
   for (int ii=0; ii<nv; ++ii)
   {
      parms->vx[ii]  = vTarget[0] - vcm.x;
      parms->vy[ii]  = vTarget[1] - vcm.y;
      parms->vz[ii]  = vTarget[2] - vcm.z;
   }

   ddcFree(vTarget);
   return parms;
}

void addVelocity(TRANSFORM* transform)
{
   ADD_VELOCITY_PARMS* parms = transform->parms;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;


   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      unsigned ind;
      if (parms->by_species)
	 ind = state->species[ii]->index;
      else
	 ind = state->group[ii]->index;
      
      state->vx[ii] += parms->vx[ind];
      state->vy[ii] += parms->vy[ind];
      state->vz[ii] += parms->vz[ind];
   }
}

THREE_VECTOR computeVcm(SYSTEM* sys)
{
   STATE* state = sys->collection->state;
   unsigned nlocal = sys->nlocal;
   
   double* vx = sys->collection->state->vx;
   double* vy = sys->collection->state->vy;
   double* vz = sys->collection->state->vz;
   
   THREE_VECTOR vmlocal; ZER03D(vmlocal);
   THREE_VECTOR vmsum;   ZER03D(vmsum);
   double mlocal = 0;
   double msum = 0;
   
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      double mass = ((ATOMTYPE_PARMS *)(state->species[ii]->parm))->mass;
      THREE_VECTOR vi;
      VSET(vi, vx[ii], vy[ii], vz[ii]);
      VSCALE(vi, mass);
      VECACUM(vmlocal, vi);
      mlocal += mass;
   }
   
   MPI_Allreduce(&vmlocal, &vmsum, 3, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   MPI_Allreduce(&mlocal, &msum, 1, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);

   VSCALE(vmsum, 1/msum);
   return vmsum;
}
