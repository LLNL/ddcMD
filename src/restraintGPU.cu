#include "restraint.h"
#include "simulate.h"
#include "nlistGPU.h"
#include "pairProcessGPU.h"
#include "gpuMemUtils.h"
#include "gpu_allocator.hpp"
#define PBC(x) x 

__global__ void restraintKernel(STATE *state, int *rback, THREE_MATRIX h0, double *e,
        int *fcx, int *fcy, int *fcz,
        double *virialxx, double *virialyy, double *virialzz,
        double * restraint_x0, double * restraint_y0, double * restraint_z0,
        double *_kb, unsigned * restraintMap, int origin, int nRestraint) 
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int r = pid;

    if (pid >= nRestraint)
    {
        return;
    }

    double *rx = state->rx;
    double *ry = state->ry;
    double *rz = state->rz;

    double *fx = state->fx;
    double *fy = state->fy;
    double *fz = state->fz;

    int ii = restraintMap[r];
    if (ii < 0)
    {
        return;
    }
    //int iii = rback[ii];

    double xBoxHalf = h0.xx / 2.0;
    double yBoxHalf = h0.yy / 2.0;
    double zBoxHalf = h0.zz / 2.0;

    double x0 = restraint_x0[r] * h0.xx;
    double y0 = restraint_y0[r] * h0.yy;
    double z0 = restraint_z0[r] * h0.zz;

    if (origin == 0)
    {
        x0 -= xBoxHalf;
        y0 -= yBoxHalf;
        z0 -= zBoxHalf;
    }

    double xDelta = rx[ii] - x0;
    double yDelta = ry[ii] - y0;
    double zDelta = rz[ii] - z0;

    if (fcx[r] > 0 && abs(xDelta) > xBoxHalf)
        PBC(preduceGPUO3(&h0, &xDelta, &yDelta, &zDelta);)

        if (fcy[r] > 0 && abs(yDelta) > yBoxHalf)
            PBC(preduceGPUO3(&h0, &xDelta, &yDelta, &zDelta);)

            if (fcz[r] > 0 && abs(zDelta) > zBoxHalf)
                PBC(preduceGPUO3(&h0, &xDelta, &yDelta, &zDelta);)

                double cxDelta = fcx[r] * xDelta;
    double cyDelta = fcy[r] * yDelta;
    double czDelta = fcz[r] * zDelta;

    double kb = _kb[r];
    double eRestraint = kb * (cxDelta * xDelta + cyDelta * yDelta + czDelta * zDelta);
    double kforce = -2 * kb;
    double fxDelta = kforce*cxDelta;
    double fyDelta = kforce*cyDelta;
    double fzDelta = kforce*czDelta;

    // printf("r %i ii %i n %i  kb %f e %f fx %f %f %f xd %f %f %f fcx %i %i %i\n", r,ii, nRestraint, kb, kb*(cxDelta*xDelta+cyDelta*yDelta+czDelta*zDelta), fxDelta, fyDelta, fzDelta, xDelta, yDelta, zDelta, fcx[r], fcy[r], fcz[r]);
    fx[ii] += fxDelta;
    fy[ii] += fyDelta;
    fz[ii] += fzDelta;
    e[ii] += eRestraint;
    virialxx[pid] -= fxDelta*cxDelta;
    virialyy[pid] -= fyDelta*cyDelta;
    virialzz[pid] -= fzDelta*czDelta;
}

void restraintGPU_parms(RESTRAINTLIST * restraint)
{
    restraint->gpu_parms = (RESTRAINTGPUPARMS*) malloc(sizeof (RESTRAINTGPUPARMS));
    RESTRAINTGPUPARMS * parms = restraint->gpu_parms;



    int nRestraint = restraint->nRestraint;
    gpu_allocator(parms->fcx, nRestraint);
    gpu_allocator(parms->fcy, nRestraint);
    gpu_allocator(parms->fcz, nRestraint);
    gpu_allocator(parms->x0, nRestraint);
    gpu_allocator(parms->y0, nRestraint);
    gpu_allocator(parms->z0, nRestraint);
    gpu_allocator(parms->kb, nRestraint);
    gpu_allocator(parms->restraintMap, nRestraint);

    int * fcx_h = (int *) malloc(sizeof (int)*nRestraint);
    int * fcy_h = (int *) malloc(sizeof (int)*nRestraint);
    int * fcz_h = (int *) malloc(sizeof (int)*nRestraint);

    double *x0 = (double *) malloc(sizeof (double)*nRestraint);
    double *y0 = (double *) malloc(sizeof (double)*nRestraint);
    double *z0 = (double *) malloc(sizeof (double)*nRestraint);
    double *kb = (double *) malloc(sizeof (double)*nRestraint);

    for (int i = 0; i < nRestraint; i++)
    {
        fcx_h[i] = restraint->restraintList[i]->fcx;
        fcy_h[i] = restraint->restraintList[i]->fcy;
        fcz_h[i] = restraint->restraintList[i]->fcz;
        x0[i] = restraint->restraintList[i]->x0;
        y0[i] = restraint->restraintList[i]->y0;
        z0[i] = restraint->restraintList[i]->z0;
        kb[i] = restraint->restraintList[i]->kb;
    }
    gpu_memcpy_host2device(parms->x0, x0, nRestraint);
    gpu_memcpy_host2device(parms->y0, y0, nRestraint);
    gpu_memcpy_host2device(parms->z0, z0, nRestraint);
    gpu_memcpy_host2device(parms->kb, kb, nRestraint);


    gpu_memcpy_host2device(parms->fcx, fcx_h, nRestraint);
    gpu_memcpy_host2device(parms->fcy, fcy_h, nRestraint);
    gpu_memcpy_host2device(parms->fcz, fcz_h, nRestraint);

}

void restraintGPU(SYSTEM *sys, RESTRAINTLIST *parms, ETYPE *e) 
{

    //printf("restraints \n");
    //STATE * state = sys->collection->state;
    SIMULATE* simulation = (SIMULATE*) sys->parent;
    DDC* ddc = simulation->ddc;
    RESTRAINTGPUPARMS * gparms = parms->gpu_parms;
    int nRestraint = parms->nRestraint;
    if (ddc->lastUpdate == sys->loop)
    {
        assignRestraintMap(sys, parms);
        gpu_memcpy_host2device(gparms->restraintMap, parms->restraintMap, nRestraint);
    }
    int blockSize = 32;
    int gridSize = ceil((float) nRestraint / blockSize);

    GPUNLIST *gnlist = sys->collection->gnlist;
    GPUVIRIALS *gpu_sion = gnlist->gpu_sion;
    restraintKernel << <gridSize, blockSize>>>(sys->collection->gpustate, gnlist->r_backbg,
        sys->box->h0, gnlist->e1,
        gparms->fcx, gparms->fcy, gparms->fcz,
        gpu_sion->xx, gpu_sion->yy, gpu_sion->zz,
        gparms->x0, gparms->y0, gparms->z0, gparms->kb,
        gparms->restraintMap, parms->origin, nRestraint);


    // SIMULATE *sim = (SIMULATE*) (sys->parent);
    //  int gpu_integrate = sim->integrator->uses_gpu;
    //  if (gpu_integrate==0) sendForceEnergyToHost(sys, e); //send energy and forces from gpu to hsot
}

