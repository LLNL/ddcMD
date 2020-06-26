#include "thermalizeTransform.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "thermalize.h"
#include "object.h"
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

typedef struct THERMALIZE_TRANSFORM_PARMS_st
{
   THERMALIZE_PARMS* parms;
}
THERMALIZE_TRANSFORM_PARMS;


static void validateObject(OBJECT* obj);


void* thermalizeTransform_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;

   char *method = NULL;

   validateObject(obj);

   THERMALIZE_TRANSFORM_PARMS* transformParms = ddcMalloc(sizeof(THERMALIZE_TRANSFORM_PARMS));

   THERMALIZE_PARMS* parms = thermalize_init();
   transformParms->parms = parms;
   
   int nGroups, nSpecies;
   species_get(NULL, NSPECIES, (void*) &nSpecies);
   group_get(NULL, NGROUPS, (void*) &nGroups);
   int ntMax = MAX(nSpecies, nGroups);
   
   double temperature[ntMax]; 
   
   int rescale;
   object_get(obj, "randomizeSeed", &parms->randomize, INT, 1, "0");
   object_get(obj, "keepVcm",       &parms->keepVcm,   INT, 1, "0");
   object_get(obj, "seed",          &parms->seed,      U64, 1, "385212586");
   object_get(obj, "rescale",       &rescale,          INT, 1, "0");
   object_get(obj, "method", &method, STRING, 1, "BOLTZMANN");

   parms->method = THERMALIZE_BOLTZMANN;

   if (strncmp(method, "BOLTZMANN", strlen("BOLTZMANN")) == 0  )
      parms->method = THERMALIZE_BOLTZMANN;

   if (strncmp(method, "RESCALE", strlen("RESCALE")) == 0  )
      parms->method = THERMALIZE_RESCALE;

   if (strncmp(method, "JUTTNER", strlen("JUTTNER")) == 0  )
      parms->method = THERMALIZE_JUTTNER;

   free(method);
   method = NULL;
   
   if (rescale != 0)
      parms->method = THERMALIZE_RESCALE;
   
   object_get(obj, "temperature", (void*)temperature, WITH_UNITS, ntMax, "0.0", "T", NULL);

   char** names = NULL;
   int nNames = 0;
   if (object_testforkeyword(obj, "species"))
   {
      nNames = object_getv(obj, "species", (void*)&names, STRING, ABORT_IF_NOT_FOUND);
      int errCode = thermalize_setTemperature(parms, THERMALIZE_BY_SPECIES, temperature, names, nNames);
      assert(errCode == 0);
   }
   else if (object_testforkeyword(obj, "groups"))
   {
      nNames = object_getv(obj, "groups", (void*)&names, STRING, ABORT_IF_NOT_FOUND);
      int errCode = thermalize_setTemperature(parms, THERMALIZE_BY_GROUP, temperature, names, nNames);
      assert(errCode == 0);
   }
   else
     thermalize_setTemperature(parms, THERMALIZE_GLOBAL, temperature, NULL, 1);

   for (int ii=0; ii<nNames; ++ii)
      ddcFree(names[ii]);
   ddcFree(names);   

   return transformParms;
}


void thermalizeTransform(TRANSFORM* transform)
{
   THERMALIZE_TRANSFORM_PARMS* parms = transform->parms;

   SIMULATE* simulate = transform->parent;

   thermalize(simulate->system, parms->parms);
}


/** Ensures that the input specification for the transform meets all
 *  requirements.  Requirements include:
 *
 *  # obsolete keyword "randomize" does not appear.
 *  # Cannot specify both group and species
 *  # If neither group nor species is specified then only a single
 *    temperature is specified
 *  # If species is specified then the number of temperatures
 *    matches the number of species
 *  # If groups is specified then the number of temperatures
 *    matches the number of groups
 *  # Any named species actually exist in the species list
 *  # Any named groups actually exist in the groups list
 */
void validateObject(OBJECT* obj)
{
   if (getRank(0) != 0)
      return;

   int nTemp = 0;
   int nsIn = 0;
   int ngIn = 0;
   double* temperature = NULL;
   char** speciesNames = NULL;
   char** groupNames   = NULL;
   nTemp = object_getv(obj, "temperature", (void*)&temperature,  WITH_UNITS, IGNORE_IF_NOT_FOUND);
   nsIn  = object_getv(obj, "species",     (void*)&speciesNames, STRING,     IGNORE_IF_NOT_FOUND);
   ngIn  = object_getv(obj, "groups",      (void*)&groupNames,   STRING,     IGNORE_IF_NOT_FOUND);
   
   if (object_testforkeyword(obj, "randomize"))
   {
      printf("FATAL: Error in thermalize transform object.\n"
	     "       You have used the obsolete keyword randomize.\n"
	     "       Use the randomizeSeed keyword instead.\n\n");
      abortAll(-1);
   }

   if (nsIn > 0 && ngIn > 0)
   {
      printf("FATAL: Error in thermalize transform object.\n"
	     "       You may not specify both groups and species\n"
	     "       choose one or the other.\n\n");
      abortAll(-1);
   }
   
   if (nsIn == 0 && ngIn == 0 && nTemp != 1)
   {
      printf("FATAL: Error in thermalize transform object.\n"
	     "       When neither groups nor species are specified\n"
	     "       you must specify exactly one temperature.\n\n");
      abortAll(-1);
   }
   
   if (nsIn > 0 && nTemp != nsIn)
   {
      printf("FATAL: Error in thermalize transform object.\n"
	     "       The number of temperatures specified does not\n"
	     "       match the number of species specified.\n\n");
      abortAll(-1);
   }
   
   if (ngIn > 0 && nTemp != ngIn)
   {
      printf("FATAL: Error in thermalize transform object.\n"
	     "       The number of temperatures specified does not\n"
	     "       match the number of groups specified.\n\n");
      abortAll(-1);
   }

   if (nsIn > 0)
   {
      for (int ii=0; ii<nsIn; ++ii)
      {
	 if (species_find(NULL, speciesNames[ii]) == NULL)
	 {
	    printf("FATAL: Error in thermalize transform object.\n"
		   "       You have specified the species name %s.\n"
		   "       There is no such species in the object database\n\n",
		   speciesNames[ii]);
	    abortAll(-1);
	 }
      }
   }
   
   if (ngIn > 0)
   {
      for (int ii=0; ii<ngIn; ++ii)
      {
	 if (group_find(NULL, groupNames[ii]) == NULL)
	 {
	    printf("FATAL: Error in thermalize transform object.\n"
		   "       You have specified the group name %s.\n"
		   "       There is no such group in the object database\n\n",
		   groupNames[ii]);
	    abortAll(-1);
	 }
      }
   }
   for (int ii=0; ii<nsIn; ++ii)
      ddcFree(speciesNames[ii]);
   for (int ii=0; ii<ngIn; ++ii)
      ddcFree(groupNames[ii]);
   ddcFree(groupNames);
   ddcFree(speciesNames);
   ddcFree(temperature);
}

