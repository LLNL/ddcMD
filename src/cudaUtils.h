#ifndef CUDAUTILS_H_
#define CUDAUTILS_H_

#if (USE_GPU==0)
#define GPROF_START(job) 
#define GPROF_END()
#define GPUFUNC(x) static inline x{}
#endif

#if (USE_GPU==1)
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <nvToolsExt.h>
//these macros allow GPU function calls
//to exist on the CPU...where they'll
//just call an empty function ie "{}"
//the inline keyword is below used
//to bypass the one declaration rule
#define GPROF_START(job) nvtxRangePushA(job)
#define GPROF_END() nvtxRangePop()
#define GPUFUNC(x) x;
#define GPUKERNEL(x)__global__  x;
#endif

#ifdef __cplusplus
extern "C" {
#endif
    GPUFUNC(void push(char * name, int c))
#ifdef __cplusplus
}
#endif

//GPU profiling functions
#if (USE_GPU==1)
#define PUSH_RANGE(name,cid) 
#define POP_RANGE() 
#else
#define PUSH_RANGE(name,cid)
#define POP_RANGE()
#endif

#define CUDA_GUARD(n) if (blockIdx.x*blockDim.x + threadIdx.x>=n) return;

#define CUDA_LOOP_X(i, n)                                 \
  for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < (n); \
       i += blockDim.x*gridDim.x) \
  

//GPU error checking functions
#define CUDA_SAFE_CALL(call) {                                          \
 cudaError_t e=call;                              \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}
#endif
