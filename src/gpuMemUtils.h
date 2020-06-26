#ifndef GPUMEMUTILS_H_
#define GPUMEMUTILS_H_
#include "cudaUtils.h"

#include "system.h"
#include "collection.h"
#include "gpunlist.h"
//these macros allow GPU function calls
//to exist on the CPU...where they'll
//just call an empty function ie "{}"
//the inline keyword is below used
//to bypass the one declaration rule
//GPUFUNC(void serializeSpecies(SYSTEM *sys, int n))
#ifdef __cplusplus
extern "C" {
#endif


    GPUFUNC(void sendGPUState(SYSTEM* sys, int n))
    GPUFUNC(void sendHostState(SYSTEM* sys, int n))
    GPUFUNC(void allocSendGPUState(COLLECTION* col, int n))
    GPUFUNC(void allocGPUnbins(GPUNLIST *gnlist, const int nBinsTot, const int blocksize))
    GPUFUNC(void allocGPUBoxInfo(SYSTEM* sys))
    GPUFUNC(void sendForceVelocityToGPU(SYSTEM *sys, int n))
    GPUFUNC(void sendForceVelocityToHost(SYSTEM *sys, int n))
    GPUFUNC(void sendPosnToHost(SYSTEM *sys, int n))
    GPUFUNC(void serializeSpecies(SYSTEM *sys, int n))
    GPUFUNC(void sendForceEnergyToHost(SYSTEM * sys, ETYPE *e))
    GPUFUNC(void deallocGPU(SYSTEM *sys, int n))

#ifdef __cplusplus
}
#endif



#endif

