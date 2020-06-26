#include "alchemyTransform.h"

#include <string.h>
#include "object.h"
#include "species.h"
#include "simulate.h"
#include "ddcMalloc.h"
#include "units.h"
#include "codata.h"
#include "box.h"

int getRank(int);


typedef struct alchemy_parms_st
{
   double volume;
   unsigned nMap;
   char** oldSpecies;
   char** newSpecies;
} ALCHEMY_PARMS;

void* alchemy_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   ALCHEMY_PARMS* parms = ddcMalloc(sizeof(ALCHEMY_PARMS));

   int nOld = object_getv(obj, "oldSpecies", (void*)&parms->oldSpecies, STRING,IGNORE_IF_NOT_FOUND);
   int nNew = object_getv(obj, "newSpecies", (void*)&parms->newSpecies, STRING,IGNORE_IF_NOT_FOUND);
   object_get(obj, "volume", &parms->volume, WITH_UNITS, 1, "0.0", "l^3", NULL);
   parms->nMap = nOld;

   if (nOld != nNew && getRank(0) == 0)
   {
      printf("Error in alchemy transform.\n"
	     "  You must specify the same number of oldSpecies and newSpecies\n");
      exit(3);
   }

   for (int ii=0; ii<nOld; ++ii)
   {
      if (species_find(NULL, parms->oldSpecies[ii]) == NULL && getRank(0) == 0)
      {
	 printf("Error in alchemy transform.\n"
		"  oldSpecies list contains species %s but there is no such\n"
		"  species in the current structure.\n", parms->oldSpecies[ii]);
	 exit(3);
      }
   }
   

   
   return parms;
}


void alchemy(TRANSFORM* transform)
{
   ALCHEMY_PARMS* parms = transform->parms;
   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;
   gid_type nGlobal = simulate->system->nglobal;
   double* rx = state->rx;
   double* ry = state->ry;
   double* rz = state->rz;

   
   SPECIES** speciesList;
   unsigned nSpecies= system_getNspecies(sys);
   species_get(NULL, SPECIESLIST, (void*)&speciesList);
   for (unsigned ii=0; ii<parms->nMap; ++ii)
   {
      for (unsigned jj=0; jj<nSpecies; ++jj)
      {
	 if (strcmp(speciesList[jj]->name, parms->oldSpecies[ii]) == 0)
	    speciesList[jj]->name = parms->newSpecies[ii];
      }
   }
   


   if (parms->volume > 0)
   {
      if (getRank(0) == 0)
      {
	 double bohrV = units_convert(parms->volume, NULL, "a0^3");
	 double angV = units_convert(parms->volume, NULL, "Ang^3");
	 printf("Rescaling simulation volume to %f bohr^3 (%f Angstrom^3) per atom\n",
		bohrV, angV);
      }
      
      double totalV = parms->volume*nGlobal;
      THREE_MATRIX hfac; 
      THREE_VECTOR rold, r; 
      box_put(sys->box, VOLUME, (void*)&totalV);
      box_get(sys->box, HFAC, (void *)&hfac);
      if ( ! matrix_equal(hfac, I_3x3))
	 for (unsigned ii=0; ii<nlocal; ++ii)
	 {
	    rold.x=rx[ii];
	    rold.y=ry[ii];
	    rold.z=rz[ii];
	    r = matrix_vector(hfac,rold); 
	    rx[ii]=r.x;
	    ry[ii]=r.y;
	    rz[ii]=r.z;
	 }
   }
   
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
