#ifndef NGLFCONSTRAINT_H
#define NGLFCONSTRAINT_H
#include "ddc.h"
#include "simulate.h"
#include "species.h"
#include "bioCharmm.h"
#include "cudaUtils.h"

typedef struct nglfconstraint_parms_st
{
    double beta, tauBarostat, T, P0;
    int isotropic;
    THREE_VECTOR *acc;
    CHARMMPOT_PARMS *charmmpot_parms;
    void (*moveFunc)(SIMULATE*simulate, STATE* state, GID_ORDER* gidOrder, int nTotal,
        char * name, RESRANGE resRange, RESI_CONN* resiConn);
    void (*volumeFunc)(COLLECTION *col, STATE *gpu_state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt, int nlocal);
    int *startsGPU;
    int *resIDGPU;
    int *starts;
    int *resID;
    int *consToResStarts;
    int *consToResStartsGPU;
    int *groupID;
    int *groupIDGPU;

} NGLFCONSTRAINT_PARMS;

NGLFCONSTRAINT_PARMS *nglfconstraint_parms(INTEGRATOR*integrator);
void nglfconstraint(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS *p);
NGLFCONSTRAINT_PARMS *nglfconstraintTest_parms(INTEGRATOR*integrator);
void nglfconstraintTest(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS *p);
NGLFCONSTRAINT_PARMS *nglfconstraintError_parms(INTEGRATOR*integrator);
void nglfconstraintError(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS *p);

#ifdef __cplusplus
extern "C"
{
#endif

GPUFUNC(void nglfconstraintGPU_parms(SYSTEM * sys, NGLFCONSTRAINT_PARMS *parms))
GPUFUNC(void nglfconstraintGPU(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS *p))
GPUFUNC(void nglfconstraintGPULangevin(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS *p))
GPUFUNC(void nglfconstraintGPULangevinLCG64(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS *p))

#ifdef __cplusplus
}
#endif



#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
