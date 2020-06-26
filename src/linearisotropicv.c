#include "linearisotropicv.h"

#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "object.h"
#include "species.h"
#include "group.h"
#include "simulate.h"
#include "ddcMalloc.h"
#include "units.h"
#include "codata.h"
#include "state.h"
#include "three_algebra.h"
#include "mpiUtils.h"
#include "random.h"


typedef struct linearisotropicv_parms_st
{
   int      by_species;
   double*  temperature;
   double*  vmax, vmin;
   LONG64   seed;
} LINEARISO_PARMS;

static void obsoleteKeywordError(const char* keyword, const char* alternate);


// i don't understand this, so i will hardcode the velocity range
// i don't understand this, so i will hardcode the velocity range
void* linearisotropicv_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   LINEARISO_PARMS* parms = ddcMalloc(sizeof(LINEARISO_PARMS));

   int nGroups, nSpecies;
   species_get(NULL, NSPECIES, (void*) &nSpecies);
   group_get(NULL, NGROUPS, (void*) &nGroups);

   int ntMax = MAX(nSpecies, nGroups);

   double* temperature = (double*) ddcMalloc(ntMax*sizeof(double));
   parms->temperature = ddcMalloc(ntMax*sizeof(double));
   for (int ii=0; ii<ntMax; ++ii)
      parms->temperature[ii] = -1.0;  // -1 -> no velocity change for this species/group.
   
   parms->by_species = -1;
   
   if (object_testforkeyword(obj, "randomize"))
      obsoleteKeywordError("randomize", "randomizeSeed");
   
   int randomize = 0;
   object_get(obj, "randomizeSeed", &randomize,   INT, 1, "0");
   object_get(obj, "seed",          &parms->seed, U64, 1, "385212586");
   
   if (randomize != 0)
      parms->seed = generateRandomSeed();
   
   char** speciesNames = NULL;
   char** groupNames = NULL;
   int nsIn = object_getv(obj, "species", (void*)&speciesNames, STRING, IGNORE_IF_NOT_FOUND);
   int ngIn = object_getv(obj, "groups",  (void*)&groupNames,   STRING, IGNORE_IF_NOT_FOUND);
   int nTemp = object_get(obj, "temperature", (void*)temperature, WITH_UNITS, ntMax, "0.0", "T", NULL);
   
   if (nsIn > 0 && ngIn > 0 )
   {
      if (getRank(0) == 0)
	 printf("Error in linearisotropicv transform object.\n"
		"Please specify either groups or species, not both.\n");
      exit(2);
   }
   
   if (nsIn == 0 && ngIn == 0 && nTemp == 1)
   {
      // handle case of only a single temperature specified with no
      // groups or species.  Make it appear that user listed all species
      // with the same temperature for each.
      nsIn = nSpecies;
      nTemp = nSpecies;
      speciesNames = (char**) ddcMalloc(nSpecies*sizeof(char*));
      for (int ii=0; ii<nSpecies; ++ii)
      {
	 temperature[ii] = temperature[0];
	 SPECIES* sPtr = species_by_index(NULL, ii);
	 speciesNames[ii] = strdup(sPtr->name);
      }
   }
   
   int tempSize = 0;
   if (nsIn > 0)
   {
      tempSize = nsIn;
      parms->by_species = 1;
   }
   if (ngIn > 0)
   {
      tempSize = ngIn;
      parms->by_species = 0;
   }
   
   if (parms->by_species < 0 )
   {
      if (getRank(0) == 0)
	 printf("Error in linearisotropicv transform object.\n"
		"Please specify either groups or species.\n");
      exit(2);
   }
   
   if (nTemp != tempSize)
   {
      if (getRank(0) == 0)
	 printf("Error in linearisotropicv transform object.\n"
		"Please provide 1 temperature for "
		"each specified group or species nTemp=%d tempSize=%d\n",nTemp,tempSize);
      exit(2);
   }
   
   for (int ii=0; ii<nsIn; ++ii)
   {
      SPECIES* species = species_find(NULL, speciesNames[ii]);
      if (species)
      {
	 int iSpecies = species->index;
	 parms->temperature[iSpecies] = temperature[ii];
      }
      ddcFree(speciesNames[ii]);
   }
   for (int ii=0; ii<ngIn; ++ii)
   {
      GROUP* group = group_find(NULL, groupNames[ii]);
      if (group)
      {
	 int iGroup = group->index;
	 parms->temperature[iGroup] = temperature[ii];
      }
      ddcFree(groupNames[ii]);
   }

   ddcFree(speciesNames);
   ddcFree(groupNames);
   ddcFree(temperature);
   
   return parms;
}


void linearisotropicv(TRANSFORM* transform)
{
   LINEARISO_PARMS* parms = transform->parms;
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
      
      double temperature = parms->temperature[ind]; // units were converted earlier.
      if (temperature < 0.0)
	 continue;  // don't change temperature of particles with negative target.
      
      PRAND48_STATE handle = prand48_init(state->label[ii], parms->seed, 0x2345612345abllu);
      double vmin = 13./.529177, vmax = 130./.529177;
      double vran = (vmax-vmin)*(prand48(&handle)) + vmin;
      
      // choose a direction, rescale velocity magnitude
      double mass = ((ATOMTYPE_PARMS *) (state->species[ii]->parm))->mass;
      double sigma = sqrt(kB*temperature/mass);
      THREE_VECTOR vv=gasdev3d0(sigma, &handle); 
      
      double dv = sqrt(vv.x*vv.x+vv.y*vv.y+vv.z*vv.z);
      vv.x*=vran/dv;
      vv.y*=vran/dv;
      vv.z*=vran/dv;
      
      state_putV(state, vv, ii); 
   }
}

static void obsoleteKeywordError(const char* keyword, const char* alternate)
{
   if (getRank(0) != 0)
      return;
   printf("FATAL: You have used the obsolete keyword %s\n"
	  "       in a linearisotropicv transform.\n", keyword);
   if (alternate != NULL)
      printf("       Use the %s keyword instead.\n", alternate);
   else
      printf("       Check the documentation for alternatives.\n");
   abortAll(-1);
}
