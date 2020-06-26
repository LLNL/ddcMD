#ifndef PAIR_PROCESS_GPU
#define PAIR_PROCESS_GPU
#include "system.h"
#include "pair.h"
#include "energyInfo.h"
#include "bioCharmmParms.h"
#include "cudaUtils.h"
//put parms in constant memory?

typedef struct ljgpu
{
    double eps;
    double sig;
    double rcut;
    double shift;
    int nspecies;
    double * charge;
} LJGPU;


#ifdef __cplusplus
extern "C"
{
#endif

void pairProcessTemplatedGpu(SYSTEM*sys, PAIR_PARMS *parms, ETYPE *e);
void pairProcessNListGpu(SYSTEM*sys, PAIR_PARMS *parms, ETYPE *e);
void binParticlesGPU(SYSTEM * sys, double rcut);
void sendForceEnergyToHost(SYSTEM *sys, ETYPE *e);
void sendEnergyToHost(SYSTEM *sys, ETYPE *e);
void permuteParticles(int *new_ids, double *rx, double *ry, double *rz, int *r_back,
                      double *rxbg, double *rybg, double *rzbg, unsigned *bins, unsigned *binsb, int n);
void permuteParticles2(COLLECTION *col, int n);
void permuteParticlesCharmm(COLLECTION *col, CHARMMPOT_PARMS *parms, int n);

GPUFUNC(void calcTensorWithGPUVirials(SYSTEM *sys, ETYPE *e))

#ifdef __cplusplus
}
#endif




#endif
