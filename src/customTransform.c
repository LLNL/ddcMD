#include "customTransform.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "ddcMalloc.h"
#include "simulate.h"
#include "object.h"
#include "codata.h"
#include "units.h"
#include "mpiUtils.h"

static void redBlueCuTransform(TRANSFORM*);
static void thermalize_delta(TRANSFORM*);
static void progressiveAlloy(TRANSFORM*);
static void shearLayerGroups(TRANSFORM*);
static void selectTransmute(TRANSFORM*);
static void grepForGid(TRANSFORM*);

/** The idea behind this file is to provide a place to create
 * "single-use" transforms---code that is to specialized to be generally
 * useful and/or isn't worth the time and trouble to do all the nice
 * input parsing etc that would be needed to make it a proper
 * transform.
 *
 * Feel free to leave whatever code you think might be helpful or useful
 * in the future in here, but please try to keep stuff in somewhat
 * organized blocks, and leave at least some comments explaining the
 * purpose and application of whatever customs you create.*/

void* customTransform_parms(TRANSFORM* transform)
{
   // This function probably stays empty because the whole idea behind
   // the custom transform is that input processing is too much of a
   // hassle.  However, if you do find yourself having to do actual work
   // in this function, please put your code in a separate function
   // (i.e., myCustom_parms), and call that function.  That will help
   // prevent confusion and conflicts.
   return NULL;
}

void customTransform(TRANSFORM* transform)
{
   // Don't add more than a function call to this function.  Delegate
   // any real work elsewhere.  Selection of which custom transform is
   // used is controlled by these hard coded if statements because it
   // prevents the compiler for issuing warnings about unused functions.
	if (0)
		redBlueCuTransform(transform);
	if (0)
		thermalize_delta(transform);
	if (0)
		progressiveAlloy(transform);
	if (0)
		shearLayerGroups(transform);
   if (0)
      selectTransmute(transform);
   if (1)
      grepForGid(transform);
}



////////////////////////////////////////////////
///    ADD YOUR CUSTOM STUFF BELOW HERE    /////
////////////////////////////////////////////////


/** Purpose written for Cu on Cu KH simuation.  Convert all atoms with
 * y>0 to CuRed and others to CuBlue.  Assumes appropriate species are
 * defined.  dfr 18-Jul-09 */
void redBlueCuTransform(TRANSFORM* transform)
{
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;

   SPECIES* CuRed = species_find(NULL, "CuRed");
   SPECIES* CuBlue = species_find(NULL, "CuBlue");
   
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      if (state->rz[ii] > 0.0)
	 state->species[ii] = CuRed;
      else
	 state->species[ii] = CuBlue;
   }
}


/** Rescale velocities to a non-Maxwellian "delta function" velocity
 *  distribution.  Each atom velocity is rescaled such that its kinetic
 *  energy is (1/2)mv^2 = (3/2)kB*T. dfr 25-Jan-11 */
void thermalize_delta(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nLocal = simulate->system->nlocal;
	double* vx = state->vx;
   double* vy = state->vy;
   double* vz = state->vz;
	double temperature;
   object_get(obj, "temperature", &temperature, WITH_UNITS, 1, "0.0", "T", NULL);
	for (unsigned ii=0; ii<nLocal; ++ii)
	{
      double mass = ((ATOMTYPE_PARMS *)(state->species[ii]->parm))->mass;
		THREE_VECTOR vi;
		VSET(vi,  vx[ii], vy[ii], vz[ii]);
		double scale = sqrt(3.0*kB*temperature/(mass*VSQ(vi)));
		VSCALE(vi, scale);
		state_putV(state, vi, ii);
	}
}

/** This transform turns atoms into species 0 with probability
 *  1-(rz/hzz) and species 1 otherwise.  Existence of two species is a
 *  prequisite.  The only point of this is to make a structure that will
 *  have interesting cgrid files that I can use to demonstrate visit.
 *  dfr 19-Apr-2011 */
void progressiveAlloy(TRANSFORM* transform)
{
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nLocal = simulate->system->nlocal;
	double* rz = state->rz;
	SPECIES** species = state->species;
	
	// No checking that box is orthorhombic
	THREE_VECTOR corner = box_get_corner(NULL);
	THREE_VECTOR diagonal = box_get_diagonal(NULL);
	double zMin = corner.z;
	double zLen = diagonal.z;

	SPECIES** speciesList = NULL;
	species_get(NULL, SPECIESLIST, &speciesList);
	int nSpecies = 0;
	species_get(NULL, NSPECIES, &nSpecies);
	assert (nSpecies>=2);
	
	for (unsigned ii=0; ii<nLocal; ++ii)
	{
		double zScaled = (rz[ii]-zMin)/zLen;
		if (zScaled > drand48())
			species[ii] = speciesList[0];
		else
			species[ii] = speciesList[1];
	}
}

/** Convert all atoms z<zMin to bottomGroup and all atoms z>zMax to
 * topGroup */
void shearLayerGroups(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nLocal = simulate->system->nlocal;
	double* rz = state->rz;
	GROUP** group = state->group;
	char* topName;
	char* botName;
	double zMin;
	double zMax;
   object_get(obj, "topGroup",    &topName, STRING, 1, "top");
   object_get(obj, "bottomGroup", &botName, STRING, 1, "bottom");
   object_get(obj, "zMin",        &zMin,    WITH_UNITS, 1, "0", "l", NULL);
   object_get(obj, "zMax",        &zMax,    WITH_UNITS, 1, "0", "l", NULL);

	GROUP* topGroup = group_find(NULL, topName);
	GROUP* botGroup = group_find(NULL, botName);
	for (unsigned ii=0; ii<nLocal; ++ii)
	{
		if (rz[ii] > zMax)
			group[ii] = topGroup;
		if (rz[ii] < zMin)
			group[ii] = botGroup;
	}
}
	
/* An embellishment of redBlueCuTransform */
/* Select particular atoms in a slab and transmute them from old to new  6/30/2012 */
void selectTransmute(TRANSFORM* transform)
{
   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system;
   STATE* state = sys->collection->state;
   unsigned nlocal = simulate->system->nlocal;

   // get box size along x, width of a slice through the box
   double sliceW = sqrt(  sys->box->h0.xx*sys->box->h0.xx
                        + sys->box->h0.xy*sys->box->h0.xy
                        + sys->box->h0.xz*sys->box->h0.xz 
                       )/6.;

   char oldSpecies[] = "proton";
   char newSpecies[] = "deuteron";
   char oldGroup  [] = "proton";
   char newGroup  [] = "deuteron";

   for (unsigned ii=0; ii<nlocal; ii++)
   {
      // compare particle ii with the list to be transformed
      // if matched, then change the species
      if (   species_find(NULL,oldSpecies) == state->species[ii]
          && group_find  (NULL,oldGroup)   == state->group  [ii]
         )
         {
          double oldMass = ( (ATOMTYPE_PARMS *)(state->species[ii]->parm) )->mass;

          // only select the appropriate atoms when inside some volume element 
          double xMin = -sliceW;
          double xMax = +sliceW;
          if (state->rx[ii]>xMin && state->rx[ii]<xMax)
             {
              state->species[ii] = species_find(NULL, newSpecies);
              state->group[ii] = group_find(NULL, newGroup);
              state->atomtype[ii] = state->species[ii]->index + (state->group[ii]->index << 16);
              double newMass = ( (ATOMTYPE_PARMS *)(state->species[ii]->parm) )->mass;

              state->vx[ii]*=sqrt(oldMass/newMass);
              state->vy[ii]*=sqrt(oldMass/newMass);
              state->vz[ii]*=sqrt(oldMass/newMass);
            }
         }
   }
}

/* Used to grab reference particle coordinate, needed to compute
   shift in shockTransform simulations */
void grepForGid(TRANSFORM* transform)
{
   const int rank = getRank(0);
   const double lc = units_convert(1.0,NULL,"Angstrom");
   OBJECT* obj = (OBJECT*) transform;
   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system;
   STATE* state = sys->collection->state;

   gid_type* gidList = NULL;
   int nGIDs = object_getv(obj, "gid", (void*)&gidList, MPI_GID_TYPE, IGNORE_IF_NOT_FOUND);
   
   if (rank == 0)
      printf("%d gids read from object.data.  Length conversion = %0.8f\n",nGIDs,lc);

   if (nGIDs > 0)
   {
      double gidZvals[nGIDs];
      for (int gg=0; gg<nGIDs; ++gg)
      {
         gidZvals[gg] = 0.0;
         for (int j=0;j<state->nlocal;j++)      
         {
            if ( state->label[j] == gidList[gg]) 
            {
               gidZvals[gg] = lc*(state->rz[j]);   // convert to Angstroms
               break ;
            }
         }
      }
      double gidZGlobal[nGIDs];
      MPI_Allreduce(&gidZvals, &gidZGlobal, nGIDs, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);

      if (rank == 0)
      {
         FILE* gidfile = fopen("gidZvals.txt","w");
         for (int gg=0; gg<nGIDs; ++gg)
            fprintf(gidfile," %"PRIu64"   %0.10f\n",gidList[gg],gidZGlobal[gg]);
         fclose(gidfile);
      }
   }
   transform->writeNeeded = 0;  // don't rewrite this structure back out
   ddcFree(gidList);
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
