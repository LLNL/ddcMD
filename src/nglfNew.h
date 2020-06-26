#ifndef NGLFNEW_H
#define NGLFNEW_H
#include "ddc.h"
#include "simulate.h"
#include "species.h"
#include "cudaUtils.h"
#include "constraint.h"

typedef struct nglfNew_parms_st 
{
    double beta,tauBarostat,P0; 
    int nConstraintGroup; 
    CONSTRAINTNEW **constraints; 

 
} NGLFNEW_PARMS;

NGLFNEW_PARMS *nglfNew_parms(INTEGRATOR*integrator);
void nglfNew(DDC *ddc, SIMULATE *simulate, NGLFNEW_PARMS *parms);

#ifdef __cplusplus
extern "C" { 
#endif

//GPUFUNC(void nglfNewGPU_parms(SYSTEM * sys, NGLFCONSTRAINT_PARMS *parms))
//GPUFUNC(void nglfconstraintGPU(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS *p))

#ifdef __cplusplus
}
#endif



#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
