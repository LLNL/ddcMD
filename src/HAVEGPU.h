#ifndef HAVEGPU_H
#include <stdio.h>
#if (USE_GPU==1)  
#define  GPUCODE(x)    x     
#include "nlistGPU.h"
#include "gpuMemUtils.h"
#include "cudaUtils.h"
#include "nglfGPU.h"
#include "bioMartiniGPU.h"
#else 
#define  GPUCODE(x)       
#endif
#endif

