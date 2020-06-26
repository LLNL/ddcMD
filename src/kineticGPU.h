#ifndef KINETIC_GPU_H
#define KINETIC_GPU_H
#include "cudaUtils.h"
#include "system.h"

#ifdef __cplusplus
extern "C" {
#endif

    GPUFUNC(void kineticGPU(SYSTEM *sys, int flag))

#ifdef __cplusplus
}
#endif


#endif
