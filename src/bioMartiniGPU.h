#ifndef BIOMARTINIGPU_H
#define BIOMARTINIGPU_H

#ifdef __cplusplus
extern "C"
{
#endif
void martiniBondGPUParms(CHARMMPOT_PARMS *parms);
void martiniNonBondGPUParms(CHARMMPOT_PARMS *parms);
void martiniGPU1(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e);
#ifdef __cplusplus
}
#endif

#endif 
