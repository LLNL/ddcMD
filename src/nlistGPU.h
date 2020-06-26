#ifndef NLISTGPU_H_
#define NLISTGPU_H_
#include "collection.h"
#include "pair.h"
#include "bioCharmmParms.h"
#include "ddc.h"
#include "cudaUtils.h"

#define PAGE_SIZE 256

//#define SHARED 128
#ifdef __cplusplus

static __device__ inline void backInBoxGPUO3(THREE_MATRIX *h1, THREE_MATRIX *h1i, double *x, double *y, double *z)
{

#define XX   0
#define YX   1
#define ZX   2
#define XY   3
#define YY   4
#define ZY   5
#define XZ   6
#define YZ   7
#define ZZ   8
#define RINT lrint
    double *h = (double *) h1;
    double *hi = (double *) h1i;
    //double da, db, dc, dx, dy, dz;
    double da, db, dc;
    da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
    db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
    dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
    *x += h[XX] * da + h[YX] * db + h[ZX] * dc;
    *y += h[XY] * da + h[YY] * db + h[ZY] * dc;
    *z += h[XZ] * da + h[YZ] * db + h[ZZ] * dc;
#undef RINT
#undef XX  
#undef YX   
#undef ZX   
#undef XY   
#undef YY   
#undef ZY   
#undef XZ   
#undef YZ   
#undef ZZ  
}

static __device__ inline void preduceGPUO3(THREE_MATRIX *h1, double *x, double *y, double *z)
{


#define XX   0
#define YX   1
#define ZX   2
#define XY   3
#define YY   4
#define ZY   5
#define XZ   6
#define YZ   7
#define ZZ   8

    double *h = (double *) h1;
    double hxx, hyy, hzz;

    hxx = h[XX];
    hyy = h[YY];
    hzz = h[ZZ];

    if (*x > 0.5 * hxx)
    {
        *x += -hxx;
    }
    if (*x < -0.5 * hxx)
    {
        *x += hxx;
    }

    if (*y > 0.5 * hyy)
    {
        *y += -hyy;
    }
    if (*y < -0.5 * hyy)
    {
        *y += hyy;
    }

    if (*z > 0.5 * hzz)
    {
        *z += -hzz;
    }
    if (*z < -0.5 * hzz)
    {
        *z += hzz;
    }

#undef XX  
#undef YX   
#undef ZX   
#undef XY   
#undef YY   
#undef ZY   
#undef XZ   
#undef YZ   
#undef ZZ   
}

static __device__ inline void preduceGPUO3x(double hxx, double hyy, double hzz, double &x, double &y, double &z)
{


#define XX   0
#define YX   1
#define ZX   2
#define XY   3
#define YY   4
#define ZY   5
#define XZ   6
#define YZ   7
#define ZZ   8

    //double *h = (double *) h1;
    //double hxx,hyy, hzz; 

    //hxx = h[XX]; hyy = h[YY]; hzz = h[ZZ]; 


    // x = x -hxx*rint(x/hxx);
    // y = y -hyy*rint(y/hyy);
    // z = z -hzz*rint(z/hzz);

    if (x > 0.5 * hxx)
    {
        x += -hxx;
    }
    if (x < -0.5 * hxx)
    {
        x += hxx;
    }

    if (y > 0.5 * hyy)
    {
        y += -hyy;
    }
    if (y < -0.5 * hyy)
    {
        y += hyy;
    }

    if (z > 0.5 * hzz)
    {
        z += -hzz;
    }
    if (z < -0.5 * hzz)
    {
        z += hzz;
    }




#undef XX  
#undef YX   
#undef ZX   
#undef XY   
#undef YY   
#undef ZY   
#undef XZ   
#undef YZ   
#undef ZZ   
}



extern "C"
{
#endif

GPUKERNEL(void unpermuteForEvalList5(int* rback, double* fx, double* fy, double* fz, double *results, 
          double* e, double *sxx, double *syy, double *szz, double *sxy, double*sxz, double *syz, int n, int nion))

void zeroGPUForceEnergyBuffers(SYSTEM *sys);
void allocPages(GPUNLIST *gnl, int particles, int density, double cutoff);
void constructList(SYSTEM *sys, double cutoff);
void pairProcessNListGpu(SYSTEM *sys, PAIR_PARMS *parms, ETYPE *e);
void charmmPairGPU(SYSTEM *sys, CHARMMPOT_PARMS *parms, ETYPE *e);
//void evalList(SYSTEM *sys, PAIR_PARMS *parms,ETYPE *e);
#ifdef __cplusplus
}
#endif

__global__ void debugPair(double *rx, double *ry, double *rz, int nion);

#endif //NLISTGPU_H_
