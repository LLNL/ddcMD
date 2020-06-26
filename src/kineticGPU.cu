#include "kineticGPU.h"
#include "system.h"
#include "three_algebra.h"
#include <cub/cub.cuh>
#include "cudaUtils.h"
#include "gpu_allocator.hpp"

typedef struct kinetic_results_str {
    double sionxx;
    double sionyy;
    double sionzz;
    double K;
} KINETIC_RESULTS;

__global__ void kineticTermsKernel(STATE *state, int *species, double *massArr, double *K,
        double *tionxx, double *tionyy, double *tionzz, int nLocal) {
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nLocal)return;
    double *vx = state->vx;
    double *vy = state->vy;
    double *vz = state->vz;

    double mass = massArr[species[pid]];
    double vxx = vx[pid] * vx[pid];
    double vyy = vy[pid] * vy[pid];
    double vzz = vz[pid] * vz[pid];
    //double vxy = vx[pid]*vy[pid]; 

    K[pid] = 0.5 * mass * (vxx + vyy + vzz);
    tionxx[pid] = mass*vxx;
    tionyy[pid] = mass*vyy;
    tionzz[pid] = mass*vzz;
    //e->tion.xy += mass*vxy;
    //e->mass += mass; 
    //e->number++; 
    //e->rk += K;
    //e->thermal_flux = J; 
    /*
         e->sion.xx -= e->tion.xx;
         e->sion.xy -= e->tion.xy;
         e->temperature = 2.0*e->rk/(3.0*nsimul); 
     */
}

void kineticGPU(SYSTEM *sys, int flag) {
    GPUNLIST *gnlist = sys->collection->gnlist;
    STATE * state = sys->collection->state;
    ETYPE * e = &sys->energyInfo;
    int nLocal = state->nlocal;

    double *tionxx = gnlist->scratch1;
    double *tionyy = gnlist->scratch2;
    double *tionzz = gnlist->scratch3;
    double *K = gnlist->scratch4;

    //calculate some temporary values that will then be reduced:
    //tionxx, tionyy, tionzz, and K (kinetic)
    int blockSize = 32;
    int gridSize = ceil((float) state->nlocal / blockSize);

    double* kinResults_d = NULL;
    gpu_allocator(kinResults_d, 4);

    size_t temp_storage_bytes = 0;
    void *d_temp_storage = NULL;
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, tionxx, kinResults_d, nLocal);
    cudaMalloc(&d_temp_storage, temp_storage_bytes);

    //reduce tionxx, tionyy, tionzz, and K
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, tionxx, kinResults_d, nLocal);
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, tionyy, kinResults_d + 1, nLocal);
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, tionzz, kinResults_d + 2, nLocal);
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, K, kinResults_d + 3, nLocal);

    KINETIC_RESULTS kr;
    KINETIC_RESULTS* krPtr = &kr;
    gpu_memcpy_device2host(krPtr, kinResults_d, 1);

    if (d_temp_storage) cudaFree(d_temp_storage);
    if (kinResults_d) cudaFree(kinResults_d);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    //ETYPE *e = &sys->energyInfo;
    //unsigned nlocal = sys->nlocal;
    gid_type nsimul = sys->nglobal;


    e->rk = kr.K;
    e->temperature = 2.0 * e->rk / (3.0 * nsimul);
}

