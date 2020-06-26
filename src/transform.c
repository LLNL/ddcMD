#include "transform.h"

#include <string.h>
#include <ctype.h>
#include "addVelocity.h"
#include "selectSubset.h"
#include "replicate.h"
#include "thermalizeTransform.h"
#include "transectMorph.h"
#include "alchemyTransform.h"
#include "boxTransform.h"
#include "linearisotropicv.h"
#include "customTransform.h"
#include "gidShuffleTransform.h"
#include "projectileTransform.h"
#include "shockTransform.h"
#include "impactTransform.h"
#include "ddcenergy.h"
#include "bisectionLoadBalance.h"

#include "object.h"
#include "simulate.h"
#include "io.h"
#include "units.h"
#include "format.h"

#include <assert.h>

int getRank(int);

TRANSFORM *transform_init(void* parent, char* name)
{
   char* type;
   char* rateString; 
   TRANSFORM* transform =
      (TRANSFORM*) object_initialize(name, "TRANSFORM", sizeof(TRANSFORM));
   transform->parent = parent;
   transform->itype = NO_TRANSFORM;
   transform->writeNeeded = 1;
   object_get((OBJECT *) transform, "type", &type, STRING, 1, "NONE");
   transform->rate = atStartThenExit; 
   object_get((OBJECT *) transform, "rate", &rateString, STRING, 1, "NONE");
   if (strcmp(rateString,"NONE") != 0) 
   {
      if (strcasecmp(rateString, "atCheckpoint") == 0) transform->rate = simulate_getCheckpointRate(NULL); 
      if (strcasecmp(rateString, "atMaxLoop") == 0)transform->rate = simulate_getMaxLoop(NULL); 
      if (isdigit(rateString[0]))
      {
         char *eptr; 
         transform->rate=strtol(rateString,&eptr, 10);
      }
   }
   
   if (strcasecmp(type, "ADDVELOCITY") == 0)
   {
      transform->itype = VELOCITY_TRANSFORM;
      transform->parms = addVelocity_parms(transform);
      transform->function = addVelocity;
   }
   if (strcasecmp(type, "SETVELOCITY") == 0)
   {
      transform->itype = VELOCITY_TRANSFORM;
      transform->parms = setVelocity_parms(transform);
      transform->function = addVelocity;
   }
   if (strcasecmp(type, "SELECTSUBSET") == 0)
   {
      transform->itype = SELECT_TRANSFORM;
      transform->parms = selectSubset_parms(transform);
      transform->function = selectSubset;
   }
   if (strcasecmp(type, "REPLICATE") == 0)
   {
      transform->itype = REPLICATE_TRANSFORM;
      transform->parms = replicate_parms(transform);
      transform->function = replicate;
   }
   if (strcasecmp(type, "THERMALIZE") == 0)
   {
      transform->itype = THERMALIZE_TRANSFORM;
      transform->parms = thermalizeTransform_parms(transform);
      transform->function = thermalizeTransform;
   }
   if (strcasecmp(type, "TRANSECTMORPH") == 0)
   {
      transform->itype = TRANSECTMORPH_TRANSFORM;
      transform->parms = transectMorph_parms(transform);
      transform->function = transectMorph;
   }
   if (strcasecmp(type, "LINEARISOTROPICV") == 0)
   {
      transform->itype = LINEARISOTROPICV_TRANSFORM;
      transform->parms = linearisotropicv_parms(transform);
      transform->function = linearisotropicv;
   }
   if (strcasecmp(type, "ALCHEMY") == 0)
   {
      transform->itype = ALCHEMY_TRANSFORM;
      transform->parms = alchemy_parms(transform);
      transform->function = alchemy;
   }
   if (strcasecmp(type, "box") == 0)
   {
      transform->itype = BOX_TRANSFORM;
      transform->parms = boxTransform_parms(transform);
      transform->function = boxTransform;
   }
   if (strcasecmp(type, "GIDSHUFFLE") == 0)
   {
      transform->itype = GID_SHUFFLE_TRANSFORM;
      transform->parms = gidShuffleTransform_parms(transform);
      transform->function = gidShuffleTransform;
   }
   if (strcasecmp(type, "PROJECTILE") == 0)
   {
      transform->itype = PROJECTILE_TRANSFORM;
      transform->parms = projectileTransform_parms(transform);
      transform->function = projectileTransform;
   }
   if (strcasecmp(type, "SHOCK") == 0)
   {
      transform->itype = SHOCK_TRANSFORM;
      transform->parms = shockTransform_parms(transform);
      transform->function = shockTransform;
   }
   if (strcasecmp(type, "APPEND") == 0)
   {
      transform->itype = APPEND_TRANSFORM;
      transform->parms = appendTransform_parms(transform);
      transform->function = appendTransform;
   }
   if (strcasecmp(type, "ASSIGNGROUPS") == 0)
   {
      transform->itype = ASSIGNGROUPS_TRANSFORM;
      transform->parms = assignGroupsTransform_parms(transform);
      transform->function = assignGroupsTransform;
   }
   if (strcasecmp(type, "IMPACT") == 0)
   {
      transform->itype = IMPACT_TRANSFORM;
      transform->parms = impactTransform_parms(transform);
      transform->function = impactTransform;
   }
   if (strcasecmp(type, "CUSTOM") == 0)
   {
      transform->itype = CUSTOM_TRANSFORM;
      transform->parms = customTransform_parms(transform);
      transform->function = customTransform;
   }

   return transform;
}

void atRateTransforms(int nTransforms, TRANSFORM** transform)
{
   SIMULATE *simulate = simulate_getSimulate(NULL); 
   SYSTEM *sys = simulate->system;
   STATE *state=sys->collection->state; 
   int callEnergy = 0; 
   for (int ii=0; ii<nTransforms; ++ii)
   {
      if ((transform[ii]->rate > 0) && TEST0(simulate->loop, transform[ii]->rate)  ) 
      {
          if (getRank(0) == 0) printf("Performing transformation %s\n",transform[ii]->name);

         transform[ii]->function(transform[ii]);
         callEnergy = 1; 
      }
   }
   if (callEnergy) 
   {
      ddc_put(simulate->ddc, DDCNLOCAL, sys->nlocal);

      if(simulate->ddc->loadBalance->itype==BISECTION) bisectionReAssign(simulate->ddc->loadBalance); 

      ddcenergy(simulate->ddc, simulate->system, 1);
      for (int kk=0; kk<sys->ngroup; kk++) 
         sys->group[kk]->Update(sys->group[kk],FRONT_TIMESTEP,state,simulate->time,0.5*simulate->dt); // Need to have the group initialize code do this step. 
   }
}
void atCheckpointTransforms(int nTransforms, TRANSFORM** transform)
{
   for (int ii=0; ii<nTransforms; ++ii)
   {
      if (transform[ii]->rate == atCheckpoint ) 
      {
         if (getRank(0) == 0) printf("Performing transformation %s\n",transform[ii]->name);
         transform[ii]->function(transform[ii]);
      }
   }
}
void atStartThenExitTransforms(int nTransforms, TRANSFORM** transform)
{
   if (nTransforms==0) return; 
   for (int ii=0; ii<nTransforms; ++ii)
   {
         if (getRank(0) == 0) printf("Performing transformation %s\n",transform[ii]->name);
      if (transform[ii]->rate == atStartThenExit ) 
      {
         if (getRank(0) == 0) printf("Performing transformation %s\n",transform[ii]->name);
         transform[ii]->function(transform[ii]);
         if (transform[ii]->writeNeeded == 0)
         {
            if (ii+1 != nTransforms)
            {
               if (getRank(0) == 0)
                  printf("The %s transform calls special code to write its "
                        "snapshot.\n  All subsequent transforms are ignored.\n",
                        transform[ii]->name);
               return;
            }
         }
      }
   }

   TRANSFORM* last = transform[nTransforms-1];
   SIMULATE* simulate = last->parent;

   if (last->writeNeeded == 0)
      return;

   // we didn't get time, loop, and writedir out of the transform
   // object during transform_init because the simulate object wasn't
   // fully initialized at the time so it would have been impossible
   // to get sensible defaut values.
   char def[32];
   OBJECT* obj = (OBJECT*) last;
   sprintf(def, "%f", simulate->time);
   object_get(obj, "time", &last->time, DOUBLE, 1, def);
   last->time = units_convert(last->time,"t",NULL); 
   sprintf(def, "%"PRId64, simulate->loop);
   object_get(obj, "loop", &last->loop, U64,    1, def);
   sprintf(def, "transform.");
   sprintf(def+strlen(def), loopFormat(), last->loop);
   object_get(obj, "writedir", &last->writedir, STRING, 1, def);

   simulate->time = last->time;
   simulate->loop = last->loop;
   CreateSnapshotdir(simulate, last->writedir);
   writeRestart(simulate,1);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
