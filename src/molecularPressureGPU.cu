#include "molecularPressureGPU.h"
#include "nlistGPU.h"
#include "simulate.h"
#include "cudaUtils.h"
#include "box.h"
#include "gpu_allocator.hpp"
#include <cub/cub.cuh>

#define CHECK1(x) 

__global__ void loadMasses(double *masses, int *species, double *speciesMasses, int nion)
{

    CUDA_LOOP_X(pid, nion)
    {
        masses[pid] = speciesMasses[species[pid]];
    }
}

__global__ void calcIntraMolecularVirial(double *rx, double *ry, double *rz, int *species,
                                         double *rxM, double *ryM, double *rzM,
                                         double *dx, double *dy, double *dz,
                                         int *atomToMoleculeID, double *atomMasses, double *moleculeMasses,
                                         int * molRefAtom, THREE_MATRIX *h1, int n)
{

    CUDA_LOOP_X(pid, n)
    {
        int mid = atomToMoleculeID[pid]; //get molecule id for this atom
        int rid = molRefAtom[mid]; //get the atom id of the molecule's 'reference atom'
        int sid = species[pid]; //species of atom

        //grab mass of atom and it's encapsulating molecule
        double molMass = moleculeMasses[mid];
        double aMass = atomMasses[sid];

        //calc position of atom with respect to molecule's reference atom's position
        double x = rx[pid] - rx[rid];
        double y = ry[pid] - ry[rid];
        double z = rz[pid] - rz[rid];

        //grab box dimensions
        double hxx = h1->xx;
        double hyy = h1->yy;
        double hzz = h1->zz;

        //periodically reduce new relative atom positions
        preduceGPUO3x(hxx, hyy, hzz, x, y, z);

        //scale positions by atom's fraction of total molecule mass
        //the sum of these scaled positions within a molecule is the molecule's center of mass
        double massFrac = aMass / molMass;

        dx[pid] = x;
        dy[pid] = y;
        dz[pid] = z;

        //calc partial centers of mass
        x *= massFrac;
        y *= massFrac;
        z *= massFrac;

        //sum and commit the molecular center of masses using atomics
        atomicAdd(rxM + mid, x);
        atomicAdd(ryM + mid, y);
        atomicAdd(rzM + mid, z);
        //rxM[pid]=x;
        //ryM[pid]=y;
        //rzM[pid]=z;
    }
    return;
}

__global__ void applyVirialCorrection(int *atomToMoleculeID, int *bin_index, double * forceBuffer,
                                      double *fx, double *fy, double *fz,
                                      double *virialx, double *virialy, double *virialz,
                                      double *dx, double *dy, double *dz,
                                      double *centerMassx, double *centerMassy, double *centerMassz, int nion)
{

    CUDA_LOOP_X(pid, nion)
    {
        //int ii = bin_index[pid]; //get the index of this particle when particles sorted in bin order
        int mid = atomToMoleculeID[pid]; //get molecule id for this atom
        double fxi = forceBuffer[7 * pid + 1] + fx[pid];
        double fyi = forceBuffer[7 * pid + 2] + fy[pid];
        double fzi = forceBuffer[7 * pid + 3] + fz[pid];
        //if (mid==0) printf("pid %i r %f dx %f dx-r %f fx %f corr %f virial %f\n",pid, centerMassx[mid], dx[pid], dx[pid]-centerMassx[mid],fxi, (dx[pid]-centerMassx[mid])*fxi);
        //subtract center of mass*force on atom from atom's virial
        //virialx[pid]+=(dx[pid]-centerMassx[mid])*fxi;
        //virialy[pid]+=(dy[pid]-centerMassy[mid])*fyi;
        //virialz[pid]+=(dz[pid]-centerMassz[mid])*fzi;
        virialx[pid] = (dx[pid] - centerMassx[mid]) * fxi;
        virialy[pid] = (dy[pid] - centerMassy[mid]) * fyi;
        virialz[pid] = (dz[pid] - centerMassz[mid]) * fzi;
    }
}

THREE_MATRIX matrix_matrixG(THREE_MATRIX b, THREE_MATRIX c)
{
    THREE_MATRIX a;
    a.xx = b.xx * c.xx + b.xy * c.yx + b.xz * c.zx;
    a.xy = b.xx * c.xy + b.xy * c.yy + b.xz * c.zy;
    a.xz = b.xx * c.xz + b.xy * c.yz + b.xz * c.zz;

    a.yx = b.yx * c.xx + b.yy * c.yx + b.yz * c.zx;
    a.yy = b.yx * c.xy + b.yy * c.yy + b.yz * c.zy;
    a.yz = b.yx * c.xz + b.yy * c.yz + b.yz * c.zz;

    a.zx = b.zx * c.xx + b.zy * c.yx + b.zz * c.zx;
    a.zy = b.zx * c.xy + b.zy * c.yy + b.zz * c.zy;
    a.zz = b.zx * c.xz + b.zy * c.yz + b.zz * c.zz;
    return a;
}

THREE_VECTOR matrix_vectorG(THREE_MATRIX m, THREE_VECTOR v)
{
    THREE_VECTOR u;
    u.x = m.xx * v.x + m.xy * v.y + m.xz * v.z;
    u.y = m.yx * v.x + m.yy * v.y + m.yz * v.z;
    u.z = m.zx * v.x + m.zy * v.y + m.zz * v.z;
    return u;
}

__global__ void adjustPosnGPU(STATE *state, THREE_MATRIX m, int nLocal)
{

    CUDA_LOOP_X(pid, nLocal)
    {
        double *rx = state->rx;
        double *ry = state->ry;
        double *rz = state->rz;

        double rxo = rx[pid];
        double ryo = ry[pid];
        double rzo = rz[pid];

        double rxn = m.xx * rxo + m.xy * ryo + m.xz*rzo;
        double ryn = m.yx * rxo + m.yy * ryo + m.yz*rzo;
        double rzn = m.zx * rxo + m.zy * ryo + m.zz*rzo;
        rx[pid] = rxn;
        ry[pid] = ryn;
        rz[pid] = rzn;
        // if (pid<10)printf("rxn %f %f %f %f %f %f\n",rxo, ryo, rzo, rxn, ryn, rzn);
    }
}

void adjustPosnCPU(STATE *state, BOX_STRUCT* box)
{
    THREE_MATRIX hfac;
    box_get(box, HFAC, (void *) &hfac);
    for (int kk = 0; kk < state->nlocal; kk++)
    {
        //THREE_VECTOR old={state->rx[kk],state->ry[kk],state->rz[kk]}; 
        THREE_VECTOR old;
        old.x = state->rx[kk];
        old.y = state->ry[kk];
        old.z = state->rz[kk];
        THREE_VECTOR newbox = matrix_vectorG(hfac, old);
        state->rx[kk] = newbox.x;
        state->ry[kk] = newbox.y;
        state->rz[kk] = newbox.z;

    }
}

void changeVolumeGPU(COLLECTION *collection, STATE *gpu_state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt, int nlocal) 
{
    GPUNLIST *gnlist=collection->gnlist;    

    double btt = beta * dt / tau;
    double Pxx, Pyy, Pzz;
    Pxx = Pyy = 0.5 * (pTensor->xx + pTensor->yy);
    Pzz = pTensor->zz;
    //Pxx = Pyy = Pzz=(1.0/3.0)*(pressureTensor->xx+pressureTensor->yy+pressureTensor->zz);   // Isotropic; 

    THREE_MATRIX lambda = mzero;
    lambda.xx = cbrt(1.0 + Pxx * btt);
    lambda.yy = cbrt(1.0 + Pyy * btt);
    //lambda.xx = 1.0;
    //lambda.yy = 1.0;
    lambda.zz = cbrt(1.0 + Pzz * btt);
    THREE_MATRIX h0 = box_get_h(box);
    //double time =system_getTime(NULL);
    //printf("lambda %f %f %f\n", lambda.xx, lambda.yy, lambda.zz);
    THREE_MATRIX h = matrix_matrixG(lambda, h0);
    box_put(box, HO, &h);

    THREE_MATRIX *h0Ptr = &(box->h0);
    THREE_MATRIX *hinvPtr = &(box->hinv);
    gpu_memcpy_host2device(gnlist->hmat_g, h0Ptr, 1);
    gpu_memcpy_host2device(gnlist->hmati_g, hinvPtr, 1);

    THREE_MATRIX hfac;
    box_get(box, HFAC, (void *) &hfac);

    //adjust positions of cpu particles
    int blockSize = 32;
    int gridSize = (nlocal + blockSize - 1) / blockSize;
    adjustPosnGPU << <gridSize, blockSize>>>(gpu_state, hfac, nlocal);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
}

void changeVolumeGPUisotropic(COLLECTION *collection, STATE *gpu_state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt, int nlocal) 
{
    GPUNLIST *gnlist=collection->gnlist; 
    double btt = beta * dt / tau;
    double Pxx, Pyy, Pzz;
    //Pxx = Pyy = 0.5*(pTensor->xx+pTensor->yy);
    //Pzz =       pTensor->zz;
    Pxx = Pyy = Pzz = (1.0 / 3.0)*(pTensor->xx + pTensor->yy + pTensor->zz); // Isotropic;

    THREE_MATRIX lambda = mzero;
    lambda.xx = cbrt(1.0 + Pxx * btt);
    lambda.yy = cbrt(1.0 + Pyy * btt);
    //lambda.xx = 1.0;
    //lambda.yy = 1.0;
    lambda.zz = cbrt(1.0 + Pzz * btt);
    THREE_MATRIX h0 = box_get_h(box);
    //double time =system_getTime(NULL);
    //printf("lambda %f %f %f\n", lambda.xx, lambda.yy, lambda.zz);
    THREE_MATRIX h = matrix_matrixG(lambda, h0);
    box_put(box, HO, &h);

    THREE_MATRIX *h0Ptr = &(box->h0);
    THREE_MATRIX *hinvPtr = &(box->hinv);
    
    gpu_memcpy_host2device(gnlist->hmat_g, h0Ptr, 1);
    gpu_memcpy_host2device(gnlist->hmati_g, hinvPtr, 1);

    THREE_MATRIX hfac;
    box_get(box, HFAC, (void *) &hfac);

    //adjust positions of cpu particles
    int blockSize = 32;
    int gridSize = (nlocal + blockSize - 1) / blockSize;
    adjustPosnGPU << <gridSize, blockSize>>>(gpu_state, hfac, nlocal);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
}

void calcMolecularPressuresGPU(SYSTEM *sys) 
{
    //grap pointers
    SIMULATE *sim = (SIMULATE*) (sys->parent);
    GPUNLIST *gnlist = sys->collection->gnlist;
    STATE * gsh = sys->collection->gpustate_h;
    //STATE * state = sys->collection->state;

    int nIon = sys->nion;
    int blockSize = 32;
    int gridSize = (nIon + blockSize - 1) / blockSize;

    //calculate molecular masses, (TODO move this to happen every redomain)
    //DDC* ddc=sim->ddc;
    if (sim->ddc->lastUpdate == sys->loop)
    { // re-count residues if redomain
        loadMasses << <gridSize, blockSize>>>(gnlist->scratch, gnlist->species_g, gnlist->mass_g, nIon);
        CUDA_SAFE_CALL(cudaPeekAtLastError());
        void *d_temp_storage = NULL;
        size_t temp_storage_bytes = 0;
        cub::DeviceSegmentedReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->scratch, gnlist->molMassesg,
                                        gnlist->nMolecules, gnlist->moleculeOffsetsg, gnlist->moleculeOffsetsg + 1);
        cudaMalloc(&d_temp_storage, temp_storage_bytes);
        cub::DeviceSegmentedReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->scratch, gnlist->molMassesg,
                                        gnlist->nMolecules, gnlist->moleculeOffsetsg, gnlist->moleculeOffsetsg + 1);
        if (d_temp_storage) cudaFree(d_temp_storage);
        CUDA_SAFE_CALL(cudaPeekAtLastError());
        cudaDeviceSynchronize();
    }

    CHECK1(
        double * tmasses = (double *) malloc(sizeof (double)* gnlist->nMolecules);
        gpu_memcpy_device2host(tmasses, gnlist->molMassesg, gnlist->nMolecules);
        for (int i = 0; i < gnlist->nMolecules; i++) {

            double mass = 0;
            int start = gnlist->moleculeOffsets[i];
            
            for (int a = gnlist->moleculeOffsets[i]; a < gnlist->moleculeOffsets[i + 1]; a++) {
                ATOMTYPE_PARMS * ap = (ATOMTYPE_PARMS*) (state->species[a]->parm);
                mass += ap->mass;
            }
            unsigned gid = (state->label[start] >> 32);
            SPECIES *species = state->species[start];

            int moleculeTypeIndex = sys->moleculeClass->speciesIndexToMoleculeIndex[species->index];
            MOLECULETYPE *type = sys->moleculeClass->moleculeTypes[moleculeTypeIndex];
            printf("name %s start %i nS %i nM %i mass %f gmas %f\n", type->name, gnlist->moleculeOffsets[i], type->nSpecies, type->nMembers, mass, tmasses[i]);

        }
    )

    cudaMemset(gnlist->partialCOMx, 0, sizeof (double)*sys->nlocal);
    cudaMemset(gnlist->partialCOMy, 0, sizeof (double)*sys->nlocal);
    cudaMemset(gnlist->partialCOMz, 0, sizeof (double)*sys->nlocal);
    //calculate partial centers of mass for each molecule
    calcIntraMolecularVirial << <gridSize, blockSize>>>(gsh->rx, gsh->ry, gsh->rz, gnlist->species_g,
        gnlist->partialCOMx, gnlist->partialCOMy, gnlist->partialCOMz,
        gnlist->scratch, gnlist->scratch1, gnlist->scratch2, gnlist->moleculeListg,
        gnlist->mass_g, gnlist->molMassesg, gnlist->molReferenceAtomsg, gnlist->hmat_g, nIon);

    CUDA_SAFE_CALL(cudaPeekAtLastError());

    applyVirialCorrection << <gridSize, blockSize>>>(gnlist->moleculeListg, gnlist->r_backg, gnlist->results, gsh->fx, gsh->fy, gsh->fz,
        gnlist->virCorx, gnlist->virCory, gnlist->virCorz,
        gnlist->scratch, gnlist->scratch1, gnlist->scratch2,
        gnlist->partialCOMx, gnlist->partialCOMy, gnlist->partialCOMz, nIon);

    CUDA_SAFE_CALL(cudaPeekAtLastError());
    return;
}

