#ifndef NGLFGPU_H_
#define NGLFGPU_H_
#include "ddc.h"
#include "simulate.h"
#include "cudaUtils.h"

typedef struct nglfgpu_parms_st
{
    double beta, tauBarostat, T, P0;
} NGLFGPU_PARMS;
#ifdef __cplusplus
extern "C"
{
#endif

GPUFUNC(void nglfGPU(DDC *ddc, SIMULATE *simulate, NGLFGPU_PARMS *p))
GPUFUNC(void nglfGPULangevin(DDC *ddc, SIMULATE *simulate, NGLFGPU_PARMS *p))

GPUFUNC(void updateStateAliasesGPU(SYSTEM *sys, unsigned *nlocal, double **rx, double **ry, double **rz, double **vx, double **vy, double **vz, double **fx, double **fy, double **fz, SPECIES ***species, gid_type **label))

GPUKERNEL(void freeVelocityUpdateGPU(int *species, double *mass, STATE *state, int nlocal, double dt))
GPUKERNEL(void freeVelocityUpdateInterleavedGPU(int *species, int *rback, double *forceBuffer, double *mass, STATE *state, int nlocal, double dt))
GPUKERNEL(void langevin_velocityUpdateGPU_FrontTimestep(int* species, STATE *state, double tau, double kBT, double *mass, int nlocal, double dt, double vcmx, double vcmy, double vcmz, double *randomx, double *randomy, double *randomz))

GPUKERNEL(void langevin_velocityUpdateGPUInterleaved_FrontTimestep(int* species, STATE *state, double *forceBuffer, double tau, double kBT, double *mass, int nlocal, double dt, double vcmx, double vcmy, double vcmz, double *randomx, double *randomy, double *randomz))

GPUKERNEL(void langevin_velocityUpdateGPU_BackTimestep(int* species, STATE *state, double tau, double kBT, double *mass, int nlocal, double dt, double vcmx, double vcmy, double vcmz, double *randomx, double *randomy, double *randomz))

GPUKERNEL(void langevin_velocityUpdateGPUInterleaved_BackTimestep(int* species, STATE *state, double *forceBuffer, double tau, double kBT, double *mass, int nlocal, double dt, double vcmx, double vcmy, double vcmz, double *randomx, double *randomy, double *randomz))

GPUKERNEL(void freePositionUpdateGPU(STATE *state, int nlocal, double dt))

#ifdef __cplusplus
}
#endif

#endif //NGLFGPU_H_
