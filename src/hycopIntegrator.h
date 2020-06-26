#ifndef HYCOPINTEGRATOR_H
#define HYCOPINTEGRATOR_H
#include "ddc.h"
#include "simulate.h"
#include "species.h"
typedef struct hycopIntegrator_parms_st
{
} HYCOPINTEGRATOR_PARMS;
HYCOPINTEGRATOR_PARMS *hycopIntegrator_parms(INTEGRATOR*integrator);
void hycopIntegrator(DDC *ddc, SIMULATE *simulate, HYCOPINTEGRATOR_PARMS *p);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
