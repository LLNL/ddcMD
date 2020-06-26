#ifndef BINPROCESS_H_
#define BINPROCESS_H_
#include "system.h"
#include "pair.h"
#ifdef __cplusplus
extern "C"
{
#endif

void testProcessor(SYSTEM* sys);
void pairProcessTemplatedLJ(SYSTEM*sys, PAIR_PARMS *parms, ETYPE *e);
#ifdef __cplusplus
}
#endif
#endif
