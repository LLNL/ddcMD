#ifndef BONDED_GPU
#define BONDED_GPU
#include "bioCharmm.h"
#include "bioCharmmParms.h"
#include "system.h"


#ifdef __cplusplus
extern "C"
{
#endif
void allocResiConManaged(SYSTEM *sys, CHARMMPOT_PARMS *parms, CHARMM_PARMS* charmmParms);
void allocResiCon(SYSTEM *sys, CHARMMPOT_PARMS *parms);
void migrateResiCon(SYSTEM *sys, CHARMMPOT_PARMS *parms);
void charmmConvalentGPU(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e);
void radix_sort(gid_type* gids, int *ids, int *buffer, int nLocal, int nIon);
#ifdef __cplusplus
}
#endif

#endif //BONDED_GPU
