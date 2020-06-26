#ifndef CUDATYPES_H_
#define CUDATYPES_H_

//#include <cuda_runtime_api.h>
//#include <cuda.h>

//#ifdef __CUDACC__
#include <cuda_fp16.h>
#include <curand.h>

//these macros allow GPU types
//to exist on the CPU...where they'll
//just substitute normal CPU types

//typedef half HALF;
//typedef curandGenerator_t GPU_RAND; 
//#endif

typedef struct gputypes_st
{
   half *rxbg_h;
   half *rybg_h;
   half *rzbg_h;
   curandGenerator_t gen;

} GPUTYPES;

#endif
