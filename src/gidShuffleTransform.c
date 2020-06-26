#include "gidShuffleTransform.h"

#include <mpi.h>

#include "transform.h"
#include "ddcMalloc.h"
#include "object.h"
#include "simulate.h"
#include "gidShuffle.h"
#include "gid.h"
#include "mpiUtils.h"

typedef struct gidShuffleTransform_parms_st
{
   int reset;
   int shuffle;
   LONG64 seed;
} GID_SHUFFLE_TRANSFORM_PARMS;

void* gidShuffleTransform_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   GID_SHUFFLE_TRANSFORM_PARMS* parms = ddcMalloc(sizeof(GID_SHUFFLE_TRANSFORM_PARMS));

   object_get(obj, "reset",   &parms->reset,   INT, 1, "0");
   object_get(obj, "shuffle", &parms->shuffle, INT, 1, "1");
   object_get(obj, "seed",    &parms->seed,    U64, 1, "0X15ED6B4FD31D1CC");
	int randomize =0;
   object_get(obj, "randomizeSeed", &randomize, INT, 1, "0");
	if (randomize != 0)
		parms->seed = generateRandomSeed();

   return parms;
}

void gidShuffleTransform(TRANSFORM* transform)
{
   GID_SHUFFLE_TRANSFORM_PARMS* parms = (GID_SHUFFLE_TRANSFORM_PARMS*) transform->parms;

   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system;
   COLLECTION* collection = sys->collection;

   if (parms->reset != 0)
      gidReset(collection);
      
   if (parms->shuffle != 0)
      gidShuffle(collection, parms->seed);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
