#include "selectSubset.h"

#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

#include "ddcMalloc.h"
#include "object.h"
#include "simulate.h"
#include "box.h"
#include "mpiUtils.h"


/** Notes:
 *
 *  In the current implementation, selectSubset always multiplies all
 *  labels by 2.  This means each select transform cancels all previous
 *  transforms.  It also means that multiple select transforms will
 *  cause the label to overflow.  We might want to find a way to control
 *  whether that shift happens.
 *
 *  There should be a way to reset the gid labels.
 *
 *  One could forsee a system where multiple selections use multiple low
 *  order bits in the label and different subsetWrites pick off
 *  different bits.
*/


enum SS_METHOD {SS_NONE, SS_ASYM_GAUSSIAN, SS_BRICK};

typedef struct subsetSelect_st
{
   enum SS_METHOD method;
   void* parms;
} SUBSET_SELECT;

typedef struct aysmGaussian_parms_st
{
   double zInt;
   double t1;
   double t2;
   double w1;
   double w2;
   double nrm1;
   double nrm2;
   LONG64 seed;
} ASYM_GAUSSIAN_PARMS;

typedef struct ss_brick_parms
{
   double xMin;
   double xMax;
   double yMin;
   double yMax;
   double zMin;
   double zMax;
} SS_BRICK_PARMS;



static void asymGaussian_parms(TRANSFORM* transform);
static unsigned asymGaussian(TRANSFORM* transform);

static void ssBrick_parms(TRANSFORM* transform);
static unsigned ssBrick(TRANSFORM* transform);


void* selectSubset_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   SUBSET_SELECT* subset = ddcMalloc(sizeof(SUBSET_SELECT));
   transform->parms = subset;
				     
   subset->method = SS_NONE;
   char* method;
   object_get(obj, "method", &method, STRING, 1, "none");

   if (strcasecmp(method, "asymgaussian") == 0)
   {
      subset->method = SS_ASYM_GAUSSIAN;
      asymGaussian_parms(transform);
   }
   if (strcasecmp(method, "brick") == 0)
   {
      subset->method = SS_BRICK;
      ssBrick_parms(transform);
   }

   if (subset->method == SS_NONE && getRank(0)==0)
   {
      printf("ERROR: The method keyword is either absent or its value is\n"
	     "       unrecognized in a subset select transform.\n"
	     "       method = %s\n", method);
      abortAll(6);
   }
   
   return subset;
}

void selectSubset(TRANSFORM* transform)
{
   SUBSET_SELECT* subset = transform->parms;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;

   for (unsigned ii=0; ii<nlocal; ++ii)
      state->label[ii] *= 2;

   unsigned nPicked=0;
   switch (subset->method)
   {
     case SS_ASYM_GAUSSIAN:
      nPicked = asymGaussian(transform);
      break;
     case SS_BRICK:
      nPicked = ssBrick(transform);
      break;
     default:
      assert(0);
   }
   unsigned long long sum;
   unsigned long long n = nPicked;
   MPI_Reduce(&n, &sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, COMM_LOCAL);
   if (getRank(0) == 0)
      printf("  %llu atoms selected\n", sum);
}

void asymGaussian_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   double pickFrac;
   double density1;
   double density2;
   double t1;
   double interface;
   object_get(obj, "frac",      &pickFrac,  DOUBLE, 1, "0.001");
   object_get(obj, "density1",  &density1,  WITH_UNITS, 1, "-1","l^-3",NULL); // in atoms/A^3
   object_get(obj, "density2",  &density2,  WITH_UNITS, 1, "-1","l^-3",NULL);
   object_get(obj, "sigma1",    &t1,        WITH_UNITS, 1, "-1","l",NULL);
   object_get(obj, "interface", &interface, DOUBLE, 1, "-1"); // fractional box units (0,1)
   if ( density1 < 0 || density2 < 0 || t1 < 0 || interface < 0 )
   {
      if (getRank(0) == 0)
	 printf("Error in TRANSFORM object.\n"
		"The asymmetric gaussian method requires that "
		"the paramters\n density1, density2, sigma1, and "
		"interface all be set.\n");
      exit(2);
   }
   
   double zSize = box_get_boundingbox(NULL).z;
   double zMax = zSize/2.0;
   double zMin = -zSize/2.0;
   double zInt = (interface - 0.5) * zSize;
   double w1 = 3.0*t1;
   double t2 = t1*density1/density2;
   double w2 = 3.0*t2;

   double spi = sqrt(atan(1.0)*4.0);
   double pickFrac1 = pickFrac * 
      ( density1 * ( zMax - zInt ) + density2 * (zInt - zMin) ) / 
      ( density1 * ( zMax - zInt ) );
   double pickFrac2 = pickFrac * 
      ( density1 * ( zMax - zInt ) + density2 * (zInt - zMin)) / 
      ( density2 * ( zInt - zMin ));
   double nrm1 = ( 1.0 / ( spi * t1 )) * pickFrac1 * ( zMax - zInt );
   double nrm2 = ( 1.0 / ( spi * t2 )) * pickFrac2 * ( zInt - zMin );

   if ( nrm1 > 1.0 || nrm2 > 1.0 ) 
   {
      if ( nrm1 > nrm2 ) 
      {
	 t1 = pickFrac1 * ( zMax - zInt ) / spi;
	 w1 = t1 * 3.0;
	 w2 = w1 * density1 / density2;
	 t2 = w2/3.0;
      }
      else 
      {
	 t2 = pickFrac2 * (zInt - zMin) / spi;
	 w2 = t2 * 3.0;
	 w1 = w2 * density2 / density1;
	 t1 = w1 / 3.0;
      }
      nrm1 = ( 1.0 / ( spi * t1 )) * pickFrac1 * ( zMax - zInt );
      nrm2 = ( 1.0 / ( spi * t2 )) * pickFrac2 * ( zInt - zMin );
   }
   LONG64 seed;
   int randomize=0;
   object_get(obj, "seed",          &seed,      U64, 1, "7254867244");
   object_get(obj, "randomizeSeed", &randomize, INT, 1, "0");
   if (randomize != 0)
      seed = generateRandomSeed();

   SUBSET_SELECT* subset = (SUBSET_SELECT*) transform->parms;
   subset->parms = ddcMalloc(sizeof(ASYM_GAUSSIAN_PARMS));
   ASYM_GAUSSIAN_PARMS* parms = (ASYM_GAUSSIAN_PARMS*) subset->parms;


   parms->zInt = zInt;
   parms->t1 = t1*t1;
   parms->t2 = t2*t2;
   parms->w1 = w1;
   parms->w2 = w2;
   parms->nrm1 = nrm1;
   parms->nrm2 = nrm2;
   parms->seed = seed;
}


unsigned asymGaussian(TRANSFORM* transform)
{
   SUBSET_SELECT* subset = transform->parms;
   ASYM_GAUSSIAN_PARMS* parms = subset->parms;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;
   
   unsigned nPicked = 0;
   double zInt = parms->zInt;
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      double z = state->rz[ii];
      if (z > zInt + parms->w1 ||
	  z < zInt - parms->w2)
	 continue;

      double dz2 = (zInt-z)*(zInt-z);
      double ff = 1;
      if (z < zInt)
	 ff = exp(-dz2/parms->t2) * parms->nrm2;
      else
	 ff = exp(-dz2/parms->t1) * parms->nrm1;
      PRAND48_STATE handle = prand48_init(state->label[ii], parms->seed, 0x93496836fed8a18dllu);
      if (ff > erand48(handle.seedVec))
      {
	 state->label[ii] |= 1;
	 ++nPicked;
      }
   }
   return nPicked;
}


/** Parameters are corners (in Angstroms) of the box to include.
 *  Defaults are (effectively) the entire domain.  */
static void ssBrick_parms(TRANSFORM* transform)
{
   SUBSET_SELECT* subset = (SUBSET_SELECT*) transform->parms;
   subset->parms = ddcMalloc(sizeof(SS_BRICK_PARMS));
   SS_BRICK_PARMS* parms = (SS_BRICK_PARMS*) subset->parms;

   OBJECT* obj = (OBJECT*) transform;
   object_get(obj, "xMin", &parms->xMin, WITH_UNITS, 1, "-1e30","l",NULL);
   object_get(obj, "xMax", &parms->xMax, WITH_UNITS, 1, "1e30","l",NULL);
   object_get(obj, "yMin", &parms->yMin, WITH_UNITS, 1, "-1e30","l",NULL);
   object_get(obj, "yMax", &parms->yMax, WITH_UNITS, 1, "1e30","l",NULL);
   object_get(obj, "zMin", &parms->zMin, WITH_UNITS, 1, "-1e30","l",NULL);
   object_get(obj, "zMax", &parms->zMax, WITH_UNITS, 1, "1e30","l",NULL);

}


unsigned ssBrick(TRANSFORM* transform)
{
   SUBSET_SELECT* subset = transform->parms;
   SS_BRICK_PARMS* parms = subset->parms;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;

   unsigned nPicked = 0;
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      double x = state->rx[ii];
      double y = state->ry[ii];
      double z = state->rz[ii];

      if ( x < parms->xMin || x > parms->xMax ||
	   y < parms->yMin || y > parms->yMax ||
	   z < parms->zMin || z > parms->zMax )
	 continue;

      state->label[ii] |= 1;
      ++nPicked;
   }
   return nPicked;
}
