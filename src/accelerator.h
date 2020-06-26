#ifndef ACCELERATOR_H
#define ACCELERATOR_H
#include "state.h"
#include "gpunlist.h"

#ifdef __cplusplus
extern "C"
{
#endif

enum ACCELERATOR_CLASS { GPU_CUDA };

typedef struct gpuCudaParms_st 
{
    int maxPagesPerParticle;
    int totbins;
    int partial_sumSize;
    int checkBounds;
    GPUNLIST *gnlist; 
    STATE *gpustate, *gpustate_h;
} GPUCUDAPARMS; 
typedef struct accelerator_st 
{
    char *name;
    char *objclass;
    char *value;
    char *type; /* model */
    void *parent;
    enum ACCELERATOR_CLASS itype;
    void *parms; 
} ACCELERATOR;

ACCELERATOR *accelerator_init(void *parent, char *name);
ACCELERATOR *accelerator_getAccelerator(ACCELERATOR *accelerator);

#ifdef __cplusplus
}
#endif

#endif /* ACCELERATOR_H */

