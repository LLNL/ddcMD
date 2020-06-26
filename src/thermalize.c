#include "thermalize.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "species.h"
#include "group.h"
#include "simulate.h"
#include "ddcMalloc.h"
#include "units.h"
#include "codata.h"
#include "state.h"
#include "three_algebra.h"
#include "random.h"
#include "mpiUtils.h"

static int getIndex(STATE* state, unsigned iAtom,
		    enum THERMALIZE_ENUM selectionMode);
static void thermalize_boltzmann(SYSTEM* system, THERMALIZE_PARMS* parms);
static void thermalize_juttner(SYSTEM* system, THERMALIZE_PARMS* parms);
static void thermalize_rescale(SYSTEM* system, THERMALIZE_PARMS* parms);

static int setTemperatureBySpecies(
   THERMALIZE_PARMS* this,
   double* temperature, char** names, unsigned nNames);
static int setTemperatureByGroup(
   THERMALIZE_PARMS* this,
   double* temperature, char** names, unsigned nNames);


THERMALIZE_PARMS* thermalize_init(void)
{
   int nGroups, nSpecies;
   species_get(NULL, NSPECIES, (void*) &nSpecies);
   group_get(NULL, NGROUPS, (void*) &nGroups);
   int ntMax = MAX(nSpecies, nGroups);


   THERMALIZE_PARMS* parms = ddcMalloc(sizeof(THERMALIZE_PARMS));
   parms->selectionMode = THERMALIZE_GLOBAL;
   parms->method = THERMALIZE_BOLTZMANN;
   parms->temperature = ddcMalloc(ntMax*sizeof(double));
   for (int ii=0; ii<ntMax; ++ii)
      // -1 -> no velocity change for this species/group.
      parms->temperature[ii] = -1.0;  
   parms->seed=0x7f9e7d6a6b2allu;
   parms->randomize = 0;
   parms->keepVcm = 0;

   return parms;
}

void thermalize_destroy(THERMALIZE_PARMS* this)
{
   ddcFree(this->temperature);
}

int thermalize_setTemperature(
   THERMALIZE_PARMS* this,
   enum THERMALIZE_ENUM mode, double* temperature,
   char** names, unsigned nNames)
{
   int errCode = 0;
   this->selectionMode = mode;
   switch (mode)
   {
     case THERMALIZE_GLOBAL:
      this->temperature[0] = temperature[0];
      break;
     case THERMALIZE_BY_SPECIES:
      errCode = setTemperatureBySpecies(this, temperature, names, nNames);
      break;
     case THERMALIZE_BY_GROUP:
      errCode = setTemperatureByGroup(this, temperature, names, nNames);
      break;
     default:
      assert(1==0);
   }
   return errCode;
}

void thermalize(SYSTEM* system, THERMALIZE_PARMS* parms)
{
   if (parms->randomize != 0)
      parms->seed = generateRandomSeed();
   switch (parms->method)
   {
     case THERMALIZE_BOLTZMANN:
      thermalize_boltzmann(system, parms);
      break;
     case THERMALIZE_RESCALE:
      thermalize_rescale(system, parms);
      break;
     case THERMALIZE_JUTTNER:
      thermalize_juttner(system, parms);
      break;
     default:
      assert(1==0);
   }
}

void thermalize_boltzmann(SYSTEM* system, THERMALIZE_PARMS* parms)
{
   STATE* state = system->collection->state;
   unsigned nlocal = system->nlocal;
   
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      int ind = getIndex(state, ii, parms->selectionMode);
      double temperature = parms->temperature[ind]; 
      // don't change temperature of particles with negative target.
      if (temperature < 0.0)
	 continue;  
      
      double mass = ((ATOMTYPE_PARMS *) (state->species[ii]->parm))->mass;
      
      LONG64 modlabel = (system->loop % 10000000) + 1;
      modlabel *= state->label[ii];
      PRAND48_STATE handle = prand48_init(modlabel, parms->seed, 0x1234512345abllu);
      double sigma = sqrt(kB*temperature/mass);
      THREE_VECTOR vv = gasdev3d0(sigma, &handle); 
      state_putV(state, vv, ii); 
   }
}

void thermalize_juttner(SYSTEM* system, THERMALIZE_PARMS* parms)
{
   STATE* state = system->collection->state;
   unsigned nlocal = system->nlocal;
   double Trel, rgamma;
   double c = hc / (2*M_PI*hbar);
   double rcsq = 1/(c*c);

   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      int ind = getIndex(state, ii, parms->selectionMode);
      double temperature = parms->temperature[ind]; 
      // don't change temperature of particles with negative target.
      if (temperature < 0.0)
	  continue;  
      
      double mass = ((ATOMTYPE_PARMS *) (state->species[ii]->parm))->mass;
      
      LONG64 modlabel = (system->loop % 10000000) + 1;
      modlabel *= state->label[ii];
      PRAND48_STATE handle = prand48_init(modlabel, parms->seed, 0x1234512345abllu);

      Trel = kB*temperature*rcsq/mass;
      THREE_VECTOR p = jutdev0(&handle,Trel);

      //jutdev0 returns gamma*beta, so get v
      rgamma = 1.0/sqrt(1.0+VSQ(p));
      VSCALE(p, c*rgamma);
      state_putV(state, p, ii);
   }
}

void thermalize_rescale(SYSTEM* system, THERMALIZE_PARMS* parms)
{
   STATE* state = system->collection->state;
   unsigned nlocal = system->nlocal;
   double* vx = state->vx;
   double* vy = state->vy;
   double* vz = state->vz;

   int nTemp;
   switch (parms->selectionMode)
   {
     case THERMALIZE_GLOBAL:
      nTemp = 1;
      break;
     case THERMALIZE_BY_SPECIES:
      species_get(NULL, NSPECIES, (void*) &nTemp);
      break;
     case THERMALIZE_BY_GROUP:
      group_get(NULL, NGROUPS, (void*) &nTemp);   
      break;
     default:
      assert(1==0);
   }

   THREE_VECTOR vcmLocal[nTemp];
   double massLocal[nTemp];
   gid_type nParticlesLocal[nTemp];
   for (int ii=0; ii<nTemp; ++ii)
   {
      VSET(vcmLocal[ii], 0.0, 0.0, 0.0);
      massLocal[ii] = 0.0;
      nParticlesLocal[ii] = 0;
   }
   
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      int ind = getIndex(state, ii, parms->selectionMode);
      assert(ind<nTemp);
      double mass = ((ATOMTYPE_PARMS *)(state->species[ii]->parm))->mass;
      THREE_VECTOR vi;
      VSET(vi, vx[ii], vy[ii], vz[ii]);
      VSCALE(vi, mass);
      VECACUM(vcmLocal[ind], vi);
      massLocal[ind] += mass;
      ++nParticlesLocal[ind];
   }

   THREE_VECTOR vcmSum[nTemp];
   double massSum[nTemp];
   gid_type nParticlesSum[nTemp];

   MPI_Allreduce(vcmLocal, vcmSum, 3*nTemp, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   MPI_Allreduce(massLocal, massSum, nTemp, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   MPI_Allreduce(nParticlesLocal, nParticlesSum, nTemp, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
   
   for (int ii=0; ii<nTemp; ++ii)
      VSCALE(vcmSum[ii], 1.0/massSum[ii]);


   double tempLocal[nTemp];
   for (int ii=0; ii<nTemp; ++ii)
      tempLocal[ii] = 0.0;

   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      unsigned ind = getIndex(state, ii, parms->selectionMode);

      double mass = ((ATOMTYPE_PARMS *)(state->species[ii]->parm))->mass;
      THREE_VECTOR vi;
      VSET(vi, vx[ii], vy[ii], vz[ii]);
      VOP1(vi, -=, vcmSum[ind]);
      tempLocal[ind] += mass * VSQ(vi);
   }

   double tempSum[nTemp];
   MPI_Allreduce(tempLocal, tempSum, nTemp, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);

   for (int ii=0; ii<nTemp; ++ii)
      tempSum[ii] /= (3.0*(nParticlesSum[ii]-1));

   double vScale[nTemp];
   THREE_VECTOR vcmNew[nTemp];
   for (int ii=0; ii<nTemp; ++ii)
   {
      if (parms->temperature[ii] > 0.0)
	 vScale[ii] = sqrt(kB*parms->temperature[ii]/tempSum[ii]);
      else
	 vScale[ii] = 1.0;
      vcmNew[ii] = vcmSum[ii];
      if (parms->keepVcm == 0)
	 VSCALE(vcmNew[ii], vScale[ii]);
//      if (getRank(0) == 0)
//      {
//	 printf("ddt: ind = %d  temperature = %f, tempSum = %f, vScale = %f vcmNew =(%f, %f, %f)\n",
//		ii, parms->temperature[ii], tempSum[ii], vScale[ii], vcmNew[ii].x,vcmNew[ii].y,vcmNew[ii].z);
//      }
   }

   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      unsigned ind = getIndex(state, ii, parms->selectionMode);

      THREE_VECTOR vi;
      VSET(vi, vx[ii], vy[ii], vz[ii]);
      VOP1(vi, -=, vcmSum[ind]);
      VSCALE(vi, vScale[ind]);
      VOP1(vi, +=, vcmNew[ind]);

      state_putV(state, vi, ii); 
   }
}

int getIndex(STATE* state, unsigned iAtom, enum THERMALIZE_ENUM selectionMode)
{
   unsigned ind;
   switch (selectionMode)
   {
     case THERMALIZE_GLOBAL:
      ind = 0;
      break;
     case THERMALIZE_BY_SPECIES:
      ind = state->species[iAtom]->index;
      break;
     case THERMALIZE_BY_GROUP:
      ind = state->group[iAtom]->index;
      break;
     default:
      assert(1==0);
   }
   return ind;
}

int setTemperatureBySpecies(
   THERMALIZE_PARMS* this,
   double* temperature, char** names, unsigned nNames)
{
   for (unsigned ii=0; ii<nNames; ++ii)
   {
      SPECIES* species = species_find(NULL, names[ii]);
      if (!species)
	 return -1;
      int iSpecies = species->index;
      this->temperature[iSpecies] = temperature[ii];
   }
   return 0;
}

int setTemperatureByGroup(
   THERMALIZE_PARMS* this,
   double* temperature, char** names, unsigned nNames)
{
   for (unsigned ii=0; ii<nNames; ++ii)
   {
      GROUP* group = group_find(NULL, names[ii]);
      if (!group)
	 return -1;
      int iGroup = group->index;
      this->temperature[iGroup] = temperature[ii];
   }
   return 0;
}
