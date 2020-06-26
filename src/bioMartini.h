#ifndef BIOMARTINI_H
#define BIOMARTINI_H

#include "species.h"
#include "bioCharmmParms.h"
#include "system.h"
//#include "cudaUtils.h"
#include "bondedGPU.h"
#include "HAVEGPU.h"

#ifdef __cplusplus
extern "C"
{
#endif

//GPUFUNC(void martiniBondGPUParms(CHARMMPOT_PARMS *parms))
//GPUFUNC(void martiniNonBondGPUParms(CHARMMPOT_PARMS *parms))
//GPUFUNC(void martiniGPU1(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e))
int getCGLJindexbySpecie(SPECIES* specie, CHARMM_PARMS *charmmParms);

#ifdef __cplusplus
}
#endif



#endif /* BIOMARTINI_H */

