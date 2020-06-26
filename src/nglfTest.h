#ifndef NGLFTEST_H
#define NGLFTEST_H
#include "ddc.h"
#include "simulate.h"
#include "species.h"
typedef struct nglfTest_parms_st
{
	THREE_VECTOR *acc; 
	unsigned subDivide; 
	double rmax; 
	double highAccuarcyDt;
	double singleDt;
} NGLFTEST_PARMS;
NGLFTEST_PARMS *nglfTest_parms(INTEGRATOR*integrator);
void nglfTest(DDC *ddc, SIMULATE *simulate, NGLFTEST_PARMS *p);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
