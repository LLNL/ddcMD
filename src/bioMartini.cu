#include "nlistGPU.h"
#include "cudaUtils.h"
#include <thrust/pair.h>
#include "pairProcessGPU.h"
#include "bioCharmmParms.h"
#include "bioCharmm.h"
#include "bioMartini.h"
#include "units.h"
#include "codata.h"
#include "gpuMemUtils.h"
#include "cudaTypes.h"
#include "simulate.h"
#include "ddcMalloc.h"
#include "system.h"
#include "gpu_allocator.hpp"
//debugging only

/*
 * Converts array of structs of LJ parms to a struct of arrays to send to gpu
 * Then sends this info to the gpu
 */
void martiniNonBondGPUParms(CHARMMPOT_PARMS *parms) 
{
    SYSTEM *sys = system_getSystem(NULL);
    GPUNLIST *gnlist=sys->collection->gnlist;
    CharmmLJGPU_PARMS* ljgpu_parms = (CharmmLJGPU_PARMS*) ddcMalloc(sizeof (CharmmLJGPU_PARMS));
    parms->gpu_ljparms_h = ljgpu_parms;

    int nspecies = parms->nspecies;
    int nspecies2 = nspecies*nspecies;
    double *rcut = (double *) ddcMalloc(sizeof (double)*nspecies2);
    double *eps = (double *) ddcMalloc(sizeof (double)*nspecies2);
    double *sigma = (double *) ddcMalloc(sizeof (double)*nspecies2);
    double *shifty = (double *) ddcMalloc(sizeof (double)*nspecies2);
    ljgpu_parms->sIndex = (int*) malloc(sizeof (int)*sys->nion);

    CharmmLJ_PARMS ** ljparms = (CharmmLJ_PARMS **) parms->parms;
    for (int i = 0; i < nspecies2; i++)
    {
        eps[i] = ljparms[i]->eps;
        sigma[i] = ljparms[i]->sigma;
        rcut[i] = ljparms[i]->rcut;
        shifty[i] = ljparms[i]->shift;
    }
    ljgpu_parms->nspecies = parms->nspecies;
    ljgpu_parms->charge = gnlist->charge_bg;
    ljgpu_parms->ke = ke;
    ljgpu_parms->crf = parms->crf;
    ljgpu_parms->krf = parms->krf;
    ljgpu_parms->iepsilon_r = 1.0 / parms->epsilon_r;

    //allocate gpu parms struct on gpu
    gpu_allocator(parms->gpu_ljparms, 1);
    //cudaMalloc((void**) &(parms->gpu_ljparms), sizeof(CharmmLJGPU_PARMS));
    int size = nspecies*nspecies;

    //allocate struct member arrays on gpu
    gpu_allocator(ljgpu_parms->shift, size);
    gpu_allocator(ljgpu_parms->eps, size);
    gpu_allocator(ljgpu_parms->sigma, size);
    gpu_allocator(ljgpu_parms->rcut, size);
    gpu_allocator(ljgpu_parms->cg_species_index, sys->nion);
    gpu_allocator(ljgpu_parms->cg_species_index_b, sys->nion);

    //memcpy struct members to gpu
    gpu_memcpy_host2device(ljgpu_parms->eps, eps, nspecies2);
    gpu_memcpy_host2device(ljgpu_parms->sigma, sigma, nspecies2);
    gpu_memcpy_host2device(ljgpu_parms->rcut, rcut, nspecies2);
    gpu_memcpy_host2device(ljgpu_parms->shift, shifty, nspecies2);

    //memcpy struct to gpu
    gpu_memcpy_host2device(parms->gpu_ljparms, ljgpu_parms, 1);
    return;
}

int cmpfunc(const void * a, const void * b) 
{

    return ( (*(MOLECULE*) a).list[0] - (*(MOLECULE*) b).list[0]);
}

void moleculeList(SYSTEM *sys) 
{
    GPUNLIST * gnlist = sys->collection->gnlist;
    //MOLECULE *molecule =  sys->moleculeClass->molecule; 
    STATE *state = sys->collection->state;
    //int nMultiMolecules =  sys->moleculeClass->nMolecules; 
    MOLECULECLASS * moleculeClass = sys->moleculeClass;
    MOLECULETYPE **moleculeTypes = moleculeClass->moleculeTypes;

    //qsort(molecule, sys->moleculeClass->nMolecules,sizeof(MOLECULE), cmpfunc);
    //for (int i =0; i <sys->nspecies; i++)printf("species %i %s \n", i, sys->species[i]->name);

    //count total molecules in system
    int nMolecules = 0;
    int cnt = 0;
    gnlist->moleculeOffsets[0] = 0;
    for (int i = 0; i < state->nlocal; i++) 
   {
        //unsigned gid= (state->label[i]>>32); 
        SPECIES *species = state->species[i];

        int moleculeTypeIndex = moleculeClass->speciesIndexToMoleculeIndex[species->index];
        MOLECULETYPE *type = moleculeTypes[moleculeTypeIndex];
        gnlist->moleculeList[i] = cnt;
        if (species == type->ownershipSpecies) {
            gnlist->molReferenceAtoms[nMolecules] = i;
            nMolecules++;
            cnt += moleculeTypes[moleculeTypeIndex]->nSpecies;
            //printf("molecule %i start %i stop %i %i\n", nMolecules-1, cnt- moleculeTypes[moleculeTypeIndex]->nSpecies, cnt, moleculeTypes[moleculeTypeIndex]->nSpecies);
            gnlist->moleculeOffsets[nMolecules] = cnt;
        }
    }
    for (int k = 0; k < nMolecules; k++) {
        int index = gnlist->moleculeOffsets[k];
        for (int i = index; i < gnlist->moleculeOffsets[k + 1]; i++) {
            gnlist->moleculeList[i] = k;
        }
    }
    gnlist->nMolecules = nMolecules;
    /*
       for (int k=0;k<nMolecules;k++)
       {
          int sindex=0;
          MOLECULETYPE  *type = molecule[k].type; 
          int a0 = type->ownershipSpeciesOffset; 
          int i0 = molecule[k].list[a0];

          gnlist->moleculeList[index] = k; 
          gnlist->molReferenceAtoms[index] = a0; 
          index++;
          sindex++;

          for (int a=0;a<type->nSpecies;a++)
             if (a != a0) gnlist->moleculeList[index++] = k;
          gnlist->moleculeOffsets[k+1] = index; 
          //printf("molecule %i %s\n", k,  molecule[k].type->name);
       }
     */
    gpu_memcpy_host2device(gnlist->moleculeListg, gnlist->moleculeList, sys->nion);
    gpu_memcpy_host2device(gnlist->moleculeOffsetsg, gnlist->moleculeOffsets, (nMolecules + 1));
    gpu_memcpy_host2device(gnlist->molReferenceAtomsg, gnlist->molReferenceAtoms, (nMolecules));
    return;
}

void martiniGPU1(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e)
{
    PUSH_RANGE("martini start", 0);
    SIMULATE *sim = (SIMULATE*) (sys->parent);
    //DDC* ddc=sim->ddc;

    int gpu_integrate = sim->integrator->uses_gpu;

    if (sim->ddc->lastUpdate == sys->loop) //re-count residues if redomain
    {
        moleculeList(sys);
    }
    if (gpu_integrate == 0) sendGPUState(sys, sys->collection->state->nion); //send particles from host to gpu
    zeroGPUForceEnergyBuffers(sys); //set force and energy buffers on gpu to zero
    PUSH_RANGE("comp martini start", 1);
    charmmPairGPU(sys, parms, e); //Calculate all nonbonded terms
    cudaStreamSynchronize(0);
    charmmConvalentGPU(sys, parms, e); //Calculate all bonded terms
    POP_RANGE();
    //cudaStreamSynchronize(0);
    if (gpu_integrate == 0) sendForceEnergyToHost(sys, e); //send energy and forces from gpu to hsot

    //sendForceEnergyToHost(sys, e); //send energy and forces from gpu to hsot
    POP_RANGE();

}

void martiniBondGPUParms(CHARMMPOT_PARMS *parms)
{
    SYSTEM *sys = NULL;
    sys = system_getSystem(sys);
    allocResiCon(sys, parms);
    migrateResiCon(sys, parms);
}

