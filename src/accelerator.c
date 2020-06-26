#include <string.h>
#include <stdlib.h>

#include "accelerator.h"
#include "object.h"
#include "ddcMalloc.h"

static ACCELERATOR *current_accelerator = NULL;

void *gpuCudaParms(OBJECT* obj)
{
   GPUCUDAPARMS *parms= ddcMalloc(sizeof(GPUCUDAPARMS)); 
    object_get( obj, "maxPagesPerParticle", &parms->maxPagesPerParticle, INT, 1, "1");
  //object_get( obj, "totbins", &parms->totbins, INT, 1, "0");
  //Force the totbins to 0 at beginning.
    parms->totbins=0;
    parms->partial_sumSize=0;
  //object_get( obj, "partial_sumSize", &parms->partial_sumSize, INT, 1, "0"); 
   return (void *)parms; 
}
ACCELERATOR *accelerator_init(void *parent, char *name) 
{
    ACCELERATOR *accelerator;
    char *type;
    
    accelerator = (ACCELERATOR *) object_initialize(name, "ACCELERATOR", sizeof (ACCELERATOR));
    current_accelerator = accelerator;
    accelerator->parent = parent;

    object_get((OBJECT *) accelerator, "type", &type, STRING, 1, "NONE");
    accelerator->type = strdup(type);
    if (strcmp(type, "CUDA") == 0) 
    {
        accelerator->itype = GPU_CUDA;
        accelerator->parms = gpuCudaParms((OBJECT*)accelerator); 
    }
    
#if 0 
    object_get((OBJECT *) accelerator, "maxPagesPerParticle", &parms->maxPagesPerParticle, INT, 1, "1");
    
    object_get((OBJECT *) accelerator, "checkBounds", &parms->checkBounds, INT, 1, "0");
    //object_get((OBJECT *) accelerator, "totbins", &parms->totbins, INT, 1, "0");
    //Force the totbins to 0 at beginning.
    accelerator->totbins = 0;
    accelerator->partial_sumSize = 0;
    //object_get((OBJECT *) accelerator, "partial_sumSize", &parms->partial_sumSize, INT, 1, "0"); 
#endif

    return accelerator;
}

ACCELERATOR *accelerator_getAccelerator(ACCELERATOR *accelerator) 
{
    if (accelerator==NULL) accelerator=current_accelerator;
    return accelerator;
}
