#include "projectileTransform.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "simulate.h"
#include "object.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"

typedef struct ProjectileTransformParms_st
{
   gid_type     projectileGid;
   char*        species;
   char*        group;
   THREE_VECTOR velocity;
   int          reset;
} PROJECTILE_TRANSFORM_PARMS;

static void validateParms(PROJECTILE_TRANSFORM_PARMS* parms);
static void projectileTransform_create(TRANSFORM* transform);
static void projectileTransform_reset(TRANSFORM* transform);


/** The projectile transform converts the species, group, and velocity
 *  of a specified particle to make it into a projectile (such as for a
 *  stopping power simulation).
 *  
 *  If reset is non-zero then the velocity of all particles of the
 *  indicated species and group (must match both) are set to the
 *  specified velocity.  This is useful if you have multiple
 *  projectiles already present in a simulation and you just want to
 *  move to a different point on the stopping curve.
 */ 
void* projectileTransform_parms(TRANSFORM* transform)
{
   PROJECTILE_TRANSFORM_PARMS* parms = ddcMalloc(sizeof(PROJECTILE_TRANSFORM_PARMS));

   OBJECT* obj = (OBJECT*) transform;

   object_get(obj, "gid",      &parms->projectileGid, U64, 1, "0");
   object_get(obj, "group",    &parms->group,    STRING, 1, "");
   object_get(obj, "reset",    &parms->reset,    INT,    1, "0");
   object_get(obj, "species",  &parms->species,  STRING, 1, "");
   object_get(obj, "velocity", &parms->velocity, WITH_UNITS, 3, "0 0 0", "l/t", NULL);
	
   validateParms(parms);
   
   return parms;
}

void projectileTransform(TRANSFORM* transform)
{
   PROJECTILE_TRANSFORM_PARMS* parms = (PROJECTILE_TRANSFORM_PARMS*) transform->parms;
	if (parms->reset == 0)
	   projectileTransform_create(transform);
	else
	   projectileTransform_reset(transform);
}

void projectileTransform_create(TRANSFORM* transform)
{
   PROJECTILE_TRANSFORM_PARMS* parms = (PROJECTILE_TRANSFORM_PARMS*) transform->parms;
   SPECIES* projectileSpecies = species_find(NULL, parms->species);
   GROUP*   projectileGroup =   group_find(NULL, parms->group);


   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system; 
   STATE* state = sys->collection->state; 
   int nLocal = sys->nlocal; 
   double* vx = state->vx; 
   double* vy = state->vy; 
   double* vz = state->vz; 
   int* atomtype = state->atomtype;
   SPECIES** species = state->species;
   GROUP**   group= state->group;
   gid_type* gid = state->label;
   int gidFound = 0;
   
   for (int ii=0; ii<nLocal; ++ii)
   {
      if (gid[ii] != parms->projectileGid)
         continue;

      gidFound = 1;
      species[ii]  = projectileSpecies;
      group[ii]    = projectileGroup;
      atomtype[ii] = species[ii]->index + (group[ii]->index << 16);
      vx[ii]       = parms->velocity.x;
      vy[ii]       = parms->velocity.y;
      vz[ii]       = parms->velocity.z;
   }
   
   int sum;
   MPI_Reduce(&gidFound, & sum, 1, MPI_INT, MPI_SUM, 0, COMM_LOCAL);
   
   if (sum == 0 && getRank(0) == 0)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       No particle with specified gid.\n\n");
      abortAll(-1);
   }
   if (sum > 1 && getRank(0) == 0)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       Multiple particles with specified gid.\n\n");
      abortAll(-1);
   }
   
}

void projectileTransform_reset(TRANSFORM* transform)
{
   PROJECTILE_TRANSFORM_PARMS* parms = (PROJECTILE_TRANSFORM_PARMS*) transform->parms;
   SPECIES* projectileSpecies = species_find(NULL, parms->species);
   GROUP*   projectileGroup =   group_find(NULL, parms->group);


   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system; 
   STATE* state = sys->collection->state; 
   int nLocal = sys->nlocal; 
   double* vx = state->vx; 
   double* vy = state->vy; 
   double* vz = state->vz; 
   SPECIES** species = state->species;
   GROUP**   group= state->group;
   int nReset= 0;
   
   for (int ii=0; ii<nLocal; ++ii)
   {
	   if (species[ii] != projectileSpecies || group[ii] != projectileGroup)
		   continue;

		++nReset;
      vx[ii] = parms->velocity.x;
      vy[ii] = parms->velocity.y;
      vz[ii] = parms->velocity.z;
   }
   
   int sum;
   MPI_Reduce(&nReset, &sum, 1, MPI_INT, MPI_SUM, 0, COMM_LOCAL);
   
   if (sum == 0 && getRank(0) == 0)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       No particles match group and species (reset option chosen).\n\n");
      abortAll(-1);
   }

	if (getRank(0) == 0)
	  printf("Reset velocity of %d projectiles\n", sum);
}

void validateParms(PROJECTILE_TRANSFORM_PARMS* parms)
{
   if (getRank(0) != 0)
      return;


   if (strlen(parms->species) == 0)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       You must specify a species\n\n.");
      abortAll(-1);
   }
   if (strlen(parms->group) == 0)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       You must specify a group\n\n.");
      abortAll(-1);
   }
   if (VSQ(parms->velocity) == 0)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       You must specify a (non-zero) velocity\n\n.");
      abortAll(-1);
   }
   
   SPECIES* projectileSpecies = species_find(NULL, parms->species);
   GROUP*   projectileGroup =   group_find(NULL, parms->group);

   if (projectileSpecies == NULL)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       There is no species named %s\n\n.",
             parms->species);
      abortAll(-1);
   }
   if (projectileGroup == NULL)
   {
      printf("FATAL: Error in projectile transform object.\n"
             "       There is no group named %s\n\n.",
             parms->group);
      abortAll(-1);
   }
   
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
