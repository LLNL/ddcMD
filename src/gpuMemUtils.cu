//#include <mpi.h>
#include <stdlib.h>
#include <nvToolsExt.h>
#include "state.h"
#include "preduce.h"
#include "gpuMemUtils.h"
#include "cudaTypes.h"
#include "cudaUtils.h"
#include <time.h>
#include "gpu_allocator.hpp"
#include "units.h"
#include "accelerator.h"
#ifdef __cplusplus
extern "C" {
#endif

void allocGPUState(COLLECTION *collection, int nIon) 
{
    GPUNLIST *gnlist = collection->gnlist = (GPUNLIST*) malloc(sizeof (GPUNLIST));
    printf("start allocs %i\n", nIon);
    STATE *state = collection->state;
    double *rxg, *ryg, *rzg, *fxg, *fyg, *fzg, *vxg, *vyg, *vzg;

    gnlist->allocFactor = 3;
    //int numAccums = 6;
    //eUDA_SAFE_CALL(cudaMalloc((void **) &gnlist->accumBuffer, 6*sizeof(double)*n);)
    // Allocate and fill host data
    //int size = n*sizeof(double);
    printf("alloc size %i\n", nIon);
    double *dummy;
    gpu_allocator(dummy, 1);
    nvtxRangePushA("init_host_data");
    //int ndebug=2;

    // Allocate device data   
    gpu_allocator(rxg, nIon);
    gpu_allocator(ryg, nIon);
    gpu_allocator(rzg, nIon);

    gpu_allocator(fxg, nIon);
    gpu_allocator(fyg, nIon);
    gpu_allocator(fzg, nIon);

    gpu_allocator(vxg, nIon);
    gpu_allocator(vyg, nIon);
    gpu_allocator(vzg, nIon);

    //allocate device data for binned particles
    gpu_allocator(gnlist->r_backg, nIon);
    gpu_allocator(gnlist->r_backbg, nIon);
    gpu_allocator(gnlist->binsg, nIon);
    gpu_allocator(gnlist->binsbg, nIon);
    gpu_allocator(gnlist->rxbg, nIon);
    gpu_allocator(gnlist->rybg, nIon);
    gpu_allocator(gnlist->rzbg, nIon);
    gpu_allocator(gnlist->e_all, nIon);
    gpu_allocator(gnlist->e1, nIon);
    gnlist->e_all_h = (double *) malloc(sizeof (double)*nIon);
    //charge
    gpu_allocator(gnlist->charge_g, nIon);
    gpu_allocator(gnlist->charge_bg, nIon);

    //species
    gpu_allocator(gnlist->species_g, nIon);
    gpu_allocator(gnlist->species_bg, nIon);
    gnlist->species_h = (int*) malloc(sizeof (int)*nIon);
    gnlist->mass_h = (double*) malloc(sizeof (double)*nIon);
    gpu_allocator(gnlist->mass_g, nIon);

    //molecules
    gnlist->moleculeList = (int *) malloc(sizeof (int)*nIon);
    gnlist->moleculeOffsets = (int *) malloc(sizeof (int)*(nIon + 1));
    gnlist->molReferenceAtoms = (int *) malloc(sizeof (int)*nIon);

    gpu_allocator(gnlist->moleculeListg, nIon);
    gpu_allocator(gnlist->moleculeOffsetsg, (nIon + 1));
    gpu_allocator(gnlist->molReferenceAtomsg, nIon);
    gpu_allocator(gnlist->molMassesg, nIon);

    //allocate scratch space for gpu
    gpu_allocator(gnlist->scratch, nIon);
    gpu_allocator(gnlist->scratch1, nIon);
    gpu_allocator(gnlist->scratch2, nIon);
    gpu_allocator(gnlist->scratch3, nIon);
    gpu_allocator(gnlist->scratch4, nIon);
    gpu_allocator(gnlist->partialCOMx, nIon);
    gpu_allocator(gnlist->partialCOMy, nIon);
    gpu_allocator(gnlist->partialCOMz, nIon);
    gpu_allocator(gnlist->virCorx, nIon);
    gpu_allocator(gnlist->virCory, nIon);
    gpu_allocator(gnlist->virCorz, nIon);


    //allocate virials structs
    //CUDA_SAFE_CALL(cudaMalloc((void **) &(collection->gpu_sion), sizeof(GPUVIRIALS)));
    gnlist->gpu_sion = (GPUVIRIALS*) malloc(sizeof (GPUVIRIALS));
    gnlist->gpu_sion_h = (GPUVIRIALS*) malloc(sizeof (GPUVIRIALS));
    gnlist->sion_h = (GPUVIRIALS*) malloc(sizeof (GPUVIRIALS));
    //virial space on gpu
    gnlist->gpu_sion_h->xx = (double *) malloc(sizeof (double)*nIon);
    gnlist->gpu_sion_h->yy = (double *) malloc(sizeof (double)*nIon);
    gnlist->gpu_sion_h->zz = (double *) malloc(sizeof (double)*nIon);
    gnlist->gpu_sion_h->xy = (double *) malloc(sizeof (double)*nIon);
    gnlist->gpu_sion_h->xz = (double *) malloc(sizeof (double)*nIon);
    gnlist->gpu_sion_h->yz = (double *) malloc(sizeof (double)*nIon);

    gnlist->sion_h->xx = (double *) malloc(sizeof (double)*nIon);
    gnlist->sion_h->yy = (double *) malloc(sizeof (double)*nIon);
    gnlist->sion_h->zz = (double *) malloc(sizeof (double)*nIon);
    gnlist->sion_h->xy = (double *) malloc(sizeof (double)*nIon);
    gnlist->sion_h->xz = (double *) malloc(sizeof (double)*nIon);
    gnlist->sion_h->yz = (double *) malloc(sizeof (double)*nIon);


    //debugging
    //size_t free_bytes, total_bytes;
    //cudaMemGetInfo(&free_bytes, &total_bytes);
    //printf("Free bytes = %f total bytes = %f\n", (double)free_bytes, (double)total_bytes);

    gpu_allocator(gnlist->gpu_sion->xx, nIon);
    gpu_allocator(gnlist->gpu_sion->yy, nIon);
    gpu_allocator(gnlist->gpu_sion->zz, nIon);
    gpu_allocator(gnlist->gpu_sion->xy, nIon);
    gpu_allocator(gnlist->gpu_sion->xz, nIon);
    gpu_allocator(gnlist->gpu_sion->yz, nIon);

    // Duplicate allocation
    //gpu_allocator(gnlist->charge_g, nIon);
    //gpu_allocator(gnlist->charge_bg, nIon);

    gpu_allocator(gnlist->minsg, 3);
    gpu_allocator(gnlist->lensg, 3);
    gpu_allocator(gnlist->nbinsg, 3);

    gpu_allocator(gnlist->nbrIds, 27);
    gpu_allocator(gnlist->nbrIdsx, 3);
    gpu_allocator(gnlist->nbrIdsy, 3);
    gpu_allocator(gnlist->nbrIdsz, 3);
    gpu_allocator(gnlist->numNbrsxyz, 3);
    gpu_allocator(gnlist->mmgpu, 6);


    //allocate random number storage 
    gnlist->gpu_types = (GPUTYPES*) malloc(sizeof (GPUTYPES));
    int curandSize = nIon;
    if (curandSize % 2 != 0) curandSize++;
    gpu_allocator(gnlist->gpuRandomx, curandSize);
    gpu_allocator(gnlist->gpuRandomy, curandSize);
    gpu_allocator(gnlist->gpuRandomz, curandSize);
    printf("CHECK: curand size nIon = %d, curandSize=%d\n", nIon, curandSize);

    //half precision allocations
    gpu_allocator(gnlist->gpu_types->rxbg_h, nIon);
    gpu_allocator(gnlist->gpu_types->rybg_h, nIon);
    gpu_allocator(gnlist->gpu_types->rzbg_h, nIon);

    gpu_allocator(gnlist->listIdsg, nIon);
    // Force null pointer at beginning for allocGPUnbins condition check
    gnlist->binCountsg = NULL;
    gnlist->binCountsg2= NULL;
    gnlist->binHeadsg  = NULL;
    gnlist->binHeadsgSave  = gnlist->binHeadsg;
    gnlist->partial_sumsg = NULL;
    nvtxRangePop();

    // Allocate helper struct on the host
    collection->gpustate_h = (STATE*) malloc(sizeof (STATE));
    STATE *gpustate_h = collection->gpustate_h;

    //populate host helper struct
    gpustate_h->rx = rxg;
    gpustate_h->ry = ryg;
    gpustate_h->rz = rzg;

    gpustate_h->fx = fxg;
    gpustate_h->fy = fyg;
    gpustate_h->fz = fzg;

    gpustate_h->vx = vxg;
    gpustate_h->vy = vyg;
    gpustate_h->vz = vzg;

    gpustate_h->nlocal = state->nlocal;
    gpustate_h->nion = state->nion;

    //Allocate device struct
    //STATE *gpustate = collection->gpustate;
    gpu_allocator(collection->gpustate, 1);
    gpu_memcpy_host2device(collection->gpustate, gpustate_h, 1);
    printf("finished allocs\n");

    //set seed for curand
    curandCreateGenerator(&gnlist->gpu_types->gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gnlist->gpu_types->gen, 1234ULL);
    //curandSetPseudoRandomGeneratorSeed(gnlist->gpu_types->gen,time(NULL)); 
}

void allocGPUnbins(GPUNLIST *gnlist, const int nBinsTot, const int blocksize)
{
    ACCELERATOR *accelerator = accelerator_getAccelerator(NULL);
    GPUCUDAPARMS * accParms = (GPUCUDAPARMS *) accelerator->parms;
    
    if(accParms->totbins<nBinsTot)
   {
        if(gnlist->binCountsg != NULL)
      {
            CUDA_SAFE_CALL(cudaFree(gnlist->binCountsg);)
        }
        if(gnlist->binCountsg2 != NULL){
            CUDA_SAFE_CALL(cudaFree(gnlist->binCountsg2);)
        }
        if(gnlist->binHeadsgSave != NULL){
            CUDA_SAFE_CALL(cudaFree(gnlist->binHeadsgSave);)
        }
        
        // allocate 2X memory for bins/partial sums
        int nBinsTot2=nBinsTot*2;
        gpu_allocator(gnlist->binCountsg, nBinsTot2);
        gpu_allocator(gnlist->binCountsg2, nBinsTot2);
        gpu_allocator(gnlist->binHeadsg, nBinsTot2);
        gnlist->binHeadsgSave = gnlist->binHeadsg;

        accParms->totbins=nBinsTot2;       

        int partial_sumSize=((nBinsTot2-1)/blocksize+1)*blocksize;

        if(accParms->partial_sumSize<partial_sumSize){
            if(gnlist->partial_sumsg != NULL){
                CUDA_SAFE_CALL(cudaFree(gnlist->partial_sumsg);)
            }       
            gpu_allocator(gnlist->partial_sumsg, partial_sumSize); 
            accParms->partial_sumSize=partial_sumSize;
        }        
        
    }
    
}

void deallocGPU(SYSTEM *sys, int n) 
{
    GPUNLIST *gnlist = sys->collection->gnlist;
    COLLECTION *collection = sys->collection;
    printf("start deallocs %i\n", n);
    //STATE *state = collection->state;   
    STATE *gpustate_h = collection->gpustate_h;


    // Allocate device data   
    CUDA_SAFE_CALL(cudaFree(gnlist->rxbg);)
    CUDA_SAFE_CALL(cudaFree(gnlist->rybg);)
    CUDA_SAFE_CALL(cudaFree(gnlist->rzbg);)

    CUDA_SAFE_CALL(cudaFree(gpustate_h->rx);)
    CUDA_SAFE_CALL(cudaFree(gpustate_h->ry);)
    CUDA_SAFE_CALL(cudaFree(gpustate_h->rz);)

    CUDA_SAFE_CALL(cudaFree(gpustate_h->fx);)
    CUDA_SAFE_CALL(cudaFree(gpustate_h->fy);)
    CUDA_SAFE_CALL(cudaFree(gpustate_h->fz);)

    CUDA_SAFE_CALL(cudaFree(gpustate_h->vx);)
    CUDA_SAFE_CALL(cudaFree(gpustate_h->vy);)
    CUDA_SAFE_CALL(cudaFree(gpustate_h->vz);)

    //allocate device data for binned particles
    //int nIon =n;
    cudaFree((gnlist->r_backg));
    cudaFree((gnlist->r_backbg));
    cudaFree((gnlist->binsg));
    cudaFree((gnlist->binsbg));
    cudaFree((gnlist->rxbg));
    cudaFree((gnlist->rybg));
    cudaFree((gnlist->rzbg));
    cudaFree((gnlist->e_all));

    //charge
    CUDA_SAFE_CALL(cudaFree((gnlist->charge_g));)
    CUDA_SAFE_CALL(cudaFree((gnlist->charge_bg));)

    //species
    cudaFree((gnlist->species_g));
    cudaFree((gnlist->species_bg));
    free(gnlist->species_h);
    free(gnlist->mass_h);
    //free(gnlist->mass_g);

    //allocate scratch space for gpu
    cudaFree((gnlist->scratch));
    cudaFree((gnlist->scratch1));
    cudaFree((gnlist->scratch2));
    cudaFree((gnlist->scratch3));
    cudaFree((gnlist->scratch4));

    //allocate virials structs
    cudaFree((gnlist->gpu_sion));
    //virial space on gpu
    cudaFree((gnlist->gpu_sion_h->xx));
    cudaFree((gnlist->gpu_sion_h->yy));
    cudaFree((gnlist->gpu_sion_h->zz));
    cudaFree((gnlist->gpu_sion_h->xy));
    cudaFree((gnlist->gpu_sion_h->xz));
    cudaFree((gnlist->gpu_sion_h->yz));


    cudaFree((gnlist->minsg));
    cudaFree((gnlist->lensg));
    cudaFree((gnlist->nbinsg));

    cudaFree((gnlist->nbrIds));
    cudaFree((gnlist->nbrIdsx));
    cudaFree((gnlist->nbrIdsy));
    cudaFree((gnlist->nbrIdsz));
    cudaFree((gnlist->numNbrsxyz));

    //FIX ME ! bin counts should be allocated later
    CUDA_SAFE_CALL(cudaFree((gnlist->binCountsg)));
    cudaFree((gnlist->binCountsg2));
    cudaFree((gnlist->binHeadsg));
    gnlist->binHeadsg=NULL;
    cudaFree((gnlist->listIdsg));
    CUDA_SAFE_CALL(cudaFree((gnlist->partial_sumsg)));

}

void allocGPUBoxInfo(SYSTEM *sys) 
{
    GPUNLIST *gnlist = sys->collection->gnlist;
    int pbc = sys->box->pbc;
    int periodicity [3] = {pbc & 1, pbc & 2, pbc & 4}; //just get a different repr. of periodic b. conditions
    printf("pbc x %i y %i z %i\n", periodicity[0], periodicity[1], periodicity[2]);
    //int *pbc_g;
    gpu_allocator(gnlist->pbc_g, 3);
    gpu_memcpy_host2device(gnlist->pbc_g, periodicity, 3);


    THREE_MATRIX *h_mat = &(sys->box->h0);
    THREE_MATRIX *h_mati = &(sys->box->hinv);

    gpu_allocator(gnlist->hmat_g, 1);
    gpu_allocator(gnlist->hmati_g, 1);
    gpu_memcpy_host2device(gnlist->hmat_g, h_mat, 1);
    gpu_memcpy_host2device(gnlist->hmati_g, h_mati, 1);

}
//

void serializeSpecies(SYSTEM *sys, int n) 
{
    GPUNLIST *gnlist = sys->collection->gnlist;
    COLLECTION *collection = sys->collection;
    //SPECIES **species = collection->state->species; 

    for (int i = 0; i < sys->nion; i++) gnlist->species_h[i] = collection->state->species[i]->index;

    for (int i = 0; i < sys->nspecies; i++)  gnlist->mass_h[i] = ((ATOMTYPE_PARMS *) (sys->species[i]->parm))->mass;

    gpu_memcpy_host2device(gnlist->mass_g, gnlist->mass_h, sys->nspecies);
    gpu_memcpy_host2device(gnlist->species_g, gnlist->species_h, n);
}

void sendForceVelocityToGPU(SYSTEM *sys, int n) 
{
    STATE * state = sys->collection->state;
    STATE *gsh = sys->collection->gpustate_h;

    gpu_memcpy_host2device(gsh->vx, state->vx, n);
    gpu_memcpy_host2device(gsh->vy, state->vy, n);
    gpu_memcpy_host2device(gsh->vz, state->vz, n);

    //CUDA_SAFE_CALL(cudaMemcpy(gsh->fx, state->fx, sizeof(double)*n, cudaMemcpyHostToDevice);) 
    //CUDA_SAFE_CALL(cudaMemcpy(gsh->fy, state->fy, sizeof(double)*n, cudaMemcpyHostToDevice);) 
    //CUDA_SAFE_CALL(cudaMemcpy(gsh->fz, state->fz, sizeof(double)*n, cudaMemcpyHostToDevice);) 
    //serializeSpecies(sys, n);
}

void sendForceVelocityToHost(SYSTEM *sys, int n) {
    STATE * state = sys->collection->state;
    STATE *gsh = sys->collection->gpustate_h;
    gpu_memcpy_device2host(state->vx, gsh->vx, n);
    gpu_memcpy_device2host(state->vy, gsh->vy, n);
    gpu_memcpy_device2host(state->vz, gsh->vz, n);


    //gpu_memcpy_device2host(state->fx,  gsh->fx, n);
    //gpu_memcpy_device2host(state->fy,  gsh->fy, n);
    //gpu_memcpy_device2host(state->fz,  gsh->fz, n);

    //CUDA_SAFE_CALL(cudaMemcpy(state->fx, gsh->fx, sizeof(double)*n, cudaMemcpyDeviceToHost);) 
    //CUDA_SAFE_CALL(cudaMemcpy(state->fy, gsh->fy, sizeof(double)*n, cudaMemcpyDeviceToHost);) 
    //CUDA_SAFE_CALL(cudaMemcpy(state->fz, gsh->fz, sizeof(double)*n, cudaMemcpyDeviceToHost);) 
}

void sendPosnToHost(SYSTEM *sys, int n) {
    STATE * state = sys->collection->state;
    STATE *gsh = sys->collection->gpustate_h;
    gpu_memcpy_device2host(state->rx, gsh->rx, n);
    gpu_memcpy_device2host(state->ry, gsh->ry, n);
    gpu_memcpy_device2host(state->rz, gsh->rz, n);

}

void sendHostState(SYSTEM *sys, int n) 
{
    COLLECTION *collection = sys->collection;
    STATE *gsh = collection->gpustate_h;
    STATE *state = collection->state;
    //int size = n*sizeof(double);

    gpu_memcpy_device2host(state->rx, gsh->rx, n);
    gpu_memcpy_device2host(state->ry, gsh->ry, n);
    gpu_memcpy_device2host(state->rz, gsh->rz, n);
}

void sendGPUState(SYSTEM *sys, int n) 
{
    // Copy arrays from host to device
    GPUNLIST *gnlist = sys->collection->gnlist;
    COLLECTION *collection = sys->collection;
    STATE *gsh = collection->gpustate_h;
    STATE *state = collection->state;
    //int size = n*sizeof(double);
    //printf("send size %i particles %i \n", size, n);
    {
        double *tmpx = new double[n], *tmpy = new double[n], *tmpz = new double[n];
        for (int i = 0; i < n; i++) {
            tmpx[i] = state->rx[i];
            tmpy[i] = state->ry[i];
            tmpz[i] = state->rz[i];
            //preduce_ortho?(Hmatrix,tmpx+i,tmpy+i,tmpz+i);
            backInBox(tmpx + i, tmpy + i, tmpz + i);
        }
        gpu_memcpy_host2device(gsh->rx, tmpx, n);
        gpu_memcpy_host2device(gsh->ry, tmpy, n);
        gpu_memcpy_host2device(gsh->rz, tmpz, n);
        delete[] tmpx;
        delete[] tmpy;
        delete[] tmpz;
    }
    gpu_memcpy_host2device(gnlist->charge_g, state->q, n);
    gsh->q = gnlist->charge_g;

    /*
         gpu_memcpy_host2device(gsh->fx,  state->fx, n);
         gpu_memcpy_host2device(gsh->fy,  state->fy, n);
         gpu_memcpy_host2device(gsh->fz,  state->fz, n);
     */
    //memcopy helper struct
    gpu_memcpy_host2device(collection->gpustate, gsh, 1);
    //CUDA_SAFE_CALL(cudaMemcpy(collection->gpu_sion, collection->gpu_sion_h, sizeof(GPUVIRIALS), cudaMemcpyHostToDevice);)
    int nIon = n;
    gpu_memcpy_host2device(gnlist->gpu_sion->xx, gnlist->gpu_sion_h->xx, nIon);
    gpu_memcpy_host2device(gnlist->gpu_sion->yy, gnlist->gpu_sion_h->yy, nIon);
    gpu_memcpy_host2device(gnlist->gpu_sion->zz, gnlist->gpu_sion_h->zz, nIon);
    gpu_memcpy_host2device(gnlist->gpu_sion->xy, gnlist->gpu_sion_h->xy, nIon);
    gpu_memcpy_host2device(gnlist->gpu_sion->xz, gnlist->gpu_sion_h->xz, nIon);
    gpu_memcpy_host2device(gnlist->gpu_sion->yz, gnlist->gpu_sion_h->yz, nIon);
    THREE_MATRIX *h_mat = &(sys->box->h0);
    gpu_memcpy_host2device(gnlist->hmat_g, h_mat, 1);

    serializeSpecies(sys, sys->nion);
    //CUDA_SAFE_CALL(cudaMemcpy(collection->species_g, collection->species_h, sizeof(int)*sys->nspecies, cudaMemcpyHostToDevice);)

    //cudaMemset(collection->e_all, 0, size);
    //cudaMemset(collection->species_g, 0, n*sizeof(int));
    //cudaMemset(collection->species_bg, 0, n*sizeof(int));

    /*
       gpu_memcpy_host2device(gsh->vx,  state->vx, n);
            gpu_memcpy_host2device(gsh->vy,  state->vy, n);
            gpu_memcpy_host2device(gsh->vz,  state->vz, n);
     */
}

void push(char * name, int cid) {
    const uint32_t colors[] = {0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff, 0x0000ffff, 0x00ff0000, 0x00ffffff};
    const int num_colors = sizeof (colors) / sizeof (uint32_t);
    int color_id = cid;
    color_id = color_id % num_colors;
    nvtxEventAttributes_t eventAttrib = {0};
    eventAttrib.version = NVTX_VERSION;
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType = NVTX_COLOR_ARGB;
    eventAttrib.color = colors[color_id];
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii = name;
    nvtxRangePushEx(&eventAttrib);
}

void allocSendGPUState(COLLECTION *collection, int n) 
{
    //we allocate 1.5 times as much memory as we estimated we need, to be safe...
    //
    //CUDA_SAFE_CALL(cudaPeekAtLastError());
    //n*=3;
    allocGPUState(collection, n);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    //sendGPUState(collection, n);  
}
#ifdef __cplusplus
}
#endif

