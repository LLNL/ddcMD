#ifndef PNGLF_H
#define PNGLF_H
#include "ddc.h"
#include "simulate.h"
#include "ewald.h"
#include "nearPair.h"
#include "nuclearRx.h"
#include "three_algebra.h"

typedef struct pnglf_parms_st
{
	POTENTIAL *potential; 
	EWALD_PARMS *ewald; 
	double smallBallRadius;
	int trFlag; 
	int analytic;
	NEARPAIR_PARMS *nearPairParms;
	NUCLEARRX *nuclearRx; 

	int nRemoteParticles;
	unsigned* rpLocalIndex;  // local indices of remote nearPair particles
	unsigned* rpRemoteIndex; // remote indices of remote nearPair particles
	int* nSend;
	int* nRecv;
	
	
} PNGLF_PARMS;

//typedef struct pnglfAnalysis_parms_str
//{
//  PNGLF_PARMS *pnglfParms;
//  NULL *nearPairList
//  unsigned nNearPair;
//  int nAccum;
//  unsigned *nNearList;
//}PNGLFANALYSIS_PARMS;

PNGLF_PARMS *pnglf_parms(INTEGRATOR*integrator);
void pnglf(DDC *ddc, SIMULATE *simulate, PNGLF_PARMS *p);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
