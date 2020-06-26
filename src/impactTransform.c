#include "impactTransform.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "ddcMalloc.h"
#include "simulate.h"
#include "object.h"
#include "codata.h"
#include "units.h"
#include "mpiUtils.h"

typedef struct ImpactTransformParms_st
{
   double       zplane;
   double       radius;
   double       r2;
   double       delta;
   double       xcenter;
   double       ycenter;
   double       zcenter;
   double       delvz;
   double       frozenz1;
   double       frozenz2;
   char*        sphereGroup;
   char*        frozenGroup;
} IMPACT_TRANSFORM_PARMS;

typedef struct shortparticle_st {unsigned domainID, type,local_index; gid_type label; double rx, ry, rz;} SHORTPARTICLE ;

void copyShortParticleToState(int n,STATE *state, SHORTPARTICLE *shortParticle);
void copyStateToShortParticle(int n,STATE *state, SHORTPARTICLE *shortParticle);
int compareShortParticlesByGid(const void*pA, const void*pB);
void particleSortinfo(char *map, int mapstride, int n);

void* impactTransform_parms(TRANSFORM* transform)
{
   IMPACT_TRANSFORM_PARMS* parms = ddcMalloc(sizeof(IMPACT_TRANSFORM_PARMS));

   OBJECT* obj = (OBJECT*) transform;

   // zplane: atoms to be removed for z > zplane
   // radius: radius of sphere (must be >= 0)
   // r2: radius of sphere squared
   // delta: shift of sphere above zplane
   // xcenter: x coordinate of sphere center
   // ycenter: y coordinate of sphere center
   // zcenter: z coordinate of sphere center
   //    zcenter = zplane + sqrt(r2) + delta
   // delvz: velocity of impacting object

   object_get(obj, "zplane",   &parms->zplane,   WITH_UNITS, 1, "0.0", "l", NULL);
   object_get(obj, "radius",   &parms->radius,   WITH_UNITS, 1, "20.0", "l", NULL);
   object_get(obj, "delta",    &parms->delta,    WITH_UNITS, 1, "10.0", "l", NULL);
   object_get(obj, "xcenter",  &parms->xcenter,  WITH_UNITS, 1, "0.0", "l", NULL);
   object_get(obj, "ycenter",  &parms->ycenter,  WITH_UNITS, 1, "0.0", "l", NULL);
   object_get(obj, "delvz",    &parms->delvz,    WITH_UNITS, 1, "0.0", "l/t", NULL);
   object_get(obj, "frozenz1", &parms->frozenz1, WITH_UNITS, 1, "0.0", "l", NULL);
   object_get(obj, "frozenz2", &parms->frozenz2, WITH_UNITS, 1, "0.0", "l", NULL);
   object_get(obj, "sphereGroup",    &parms->sphereGroup,    STRING, 1, "S3");
   object_get(obj, "frozenGroup",    &parms->frozenGroup,    STRING, 1, "FR");

   if (parms->radius < 0.0)
   {
      printf("Error in impact transform.\n"
	     "  The sphere radius (radius) must be non-negative.\n");
      exit(3);
   }

   parms->r2 = parms->radius*parms->radius;                         // radius of sphere squared
   parms->zcenter = parms->zplane + sqrt(parms->r2) + parms->delta; // z coordinate of sphere center

   return parms;
}



static int  markForDeletionImpact(TRANSFORM* transform, STATE *state, double zplane, double r2, double xx, double yy, double zz)
{
   double zrel,dx,dy,dz;
   unsigned nLocal =state->nlocal; 
   for (int j =0;j<state->nlocal;j++)
   {
      zrel = state->rz[j] - zplane;
      dx = state->rx[j] - xx;
      dy = state->ry[j] - yy;
      dz = state->rz[j] - zz; 
      if ( (dx*dx + dy*dy + dz*dz) > r2 && zrel > 0.0 )
           {
               state->label[j] = (1LL<<62)-1;
               nLocal--; 
           }
   }
   return nLocal; 
}

/** Create a impact configuration. 
 *    In particular, remove all of the atoms above a plane in the simulation
 *    except for those atoms in a specified sphere.
 *    Give the sphere a specified velocity. rer Jan. 21, 2016. */
void impactTransform(TRANSFORM* transform)
{
   IMPACT_TRANSFORM_PARMS* parms = transform->parms;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   SYSTEM *sys = system_getSystem(NULL);

// parameters  lengths in Bohr; velocity in Bohr/fs
   double zplane    = parms->zplane;
   double frozenz1  = parms->frozenz1;
   double frozenz2  = parms->frozenz2;
   double r2        = parms->r2;
   double delvz     = parms->delvz;
   double xcenter   = parms->xcenter;
   double ycenter   = parms->ycenter;
   double zcenter   = parms->zcenter;

   char*  sphereGroup = parms->sphereGroup;
   GROUP* sphereGrp = group_find(NULL, sphereGroup);
   char*  frozenGroup = parms->frozenGroup;
   GROUP* frozenGrp = group_find(NULL, frozenGroup);

   int id = getRank(0);

   if ( id == 0 ) 
    printf("TRANSFORM -- impactTransform\n");

   unsigned nlocalUpdated = markForDeletionImpact(transform,state,zplane,r2,xcenter,ycenter,zcenter);

   unsigned nSubLocal = state->nlocal-nlocalUpdated;
   SHORTPARTICLE shortParticle[state->nlocal]; 
   copyStateToShortParticle(state->nlocal,state,shortParticle); 
   qsort(shortParticle,state->nlocal,sizeof(SHORTPARTICLE),compareShortParticlesByGid);
   copyShortParticleToState(state->nlocal,state,shortParticle); 
   particleSortinfo((char *)&(shortParticle->local_index),sizeof(SHORTPARTICLE),state->nlocal);
   sys->nlocal=sys->collection->size=state->nlocal=nlocalUpdated;
   gid_type nglobal = state->nlocal; 
   assert(sizeof(nglobal) == sizeof(uint64_t));
   MPI_Allreduce(&nglobal, &sys->nglobal, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, COMM_LOCAL);

   unsigned nlocal = simulate->system->nlocal;

   // shift velocity and update group name
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      if (state->rz[ii] > zplane) 
	{
           state->vz[ii] += delvz;
           state->group[ii] = sphereGrp;
	}
      if (frozenz1 <= state->rz[ii] && state->rz[ii] < frozenz2) 
           state->group[ii] = frozenGrp;
   }

   if ( nSubLocal != 0 )
    printf("rank = %d deleted atoms=%d atoms remaining = %d\n",id,nSubLocal,nlocal);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
