#ifndef MOL_PRESSURE_GPU
#define MOL_PRESSURE_GPU
#include "system.h"
#include "pair.h"
#include "energyInfo.h"
#include "bioCharmmParms.h"

#ifdef __cplusplus
extern "C"
{
#endif

GPUFUNC(void calcMolecularPressuresGPU(SYSTEM *sys))
GPUFUNC(void changeVolumeGPU(COLLECTION * col, STATE *state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt, int nlocal))
GPUFUNC(void changeVolumeGPUisotropic(COLLECTION * col, STATE *state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt, int nlocal))

#ifdef __cplusplus
}
#endif

#endif
