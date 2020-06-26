#ifndef NGLF_H
#define NGLF_H
#include "ddc.h"
#include "simulate.h"
#include "species.h"
typedef struct nglf_parms_st
{
//    double beta,tauBarostat,T,P0;
//	THREE_VECTOR *acc; 
} NGLF_PARMS;
NGLF_PARMS *nglf_parms(INTEGRATOR*integrator);
void nglf(DDC *ddc, SIMULATE *simulate, NGLF_PARMS *p);
NGLF_PARMS *nglfTest_parms(INTEGRATOR*integrator);
void nglfTest(DDC *ddc, SIMULATE *simulate, NGLF_PARMS *p);
NGLF_PARMS *nglfError_parms(INTEGRATOR*integrator);
void nglfError(DDC *ddc, SIMULATE *simulate, NGLF_PARMS *p);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
