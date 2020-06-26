#ifndef DDCENERGY_H
#define DDCENERGY_H

#include "ddc.h"
#include "system.h"
#include "energyInfo.h"
#ifdef __cplusplus
extern "C" { 
void fsumX(SYSTEM*sys, THREE_SMATRIX sion);
#else
void fsumX(SYSTEM*sys, THREE_SMATRIX*restrict sion);
#endif

int ddcenergy(DDC*ddc, SYSTEM*sys, int e_eval_flag);
void ddcUpdateTables(DDC*ddc, SYSTEM*system);
void tempstorage(DDC*ddc, SYSTEM*sys);
int  ddcUpdateAll(DDC*ddc, SYSTEM*sys, ETYPE *e, int ForcedUpdate);
int  ddcUpdateAll_(DDC*ddc, SYSTEM*sys, ETYPE *e, int ForcedUpdate);
int  ddcUpdateAll_pair(DDC*ddc, SYSTEM*sys, ETYPE *e, int ForcedUpdate);
void zeroAll(SYSTEM*sys);
void kinetic_terms(SYSTEM*sys, int flag);
void cutoffs(DDC*ddc, SYSTEM*sys);
#endif
#ifdef __cplusplus
} 
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
