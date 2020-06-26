#include "nglfconstraint.h"
//#include "nglfconstraintGPU.h"
#include "nglf.h"
#include "bioCharmmParms.h"
#include "bioCharmm.h"
#include "bioCharmmCovalent.h"
#include "nglfGPU.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include "three_algebra.h"
#include "object.h"
#include "ddc.h"
#include "box.h"
#include "species.h"
#include "ddcMalloc.h"
#include "ddcenergy.h"
#include "expandbuffer.h"
#include "auxNeighbor.h"
#include "preduce.h"
#include "group.h"
#include "energyInfo.h"
#include "cudaUtils.h"
#include "nlistGPU.h"
#include "bioGid.h"
//#include "velocityAutocorrelation.h"
#include "gpuMemUtils.h"
#include "pairProcessGPU.h"
#include "molecularPressureGPU.h"
#include "cudaTypes.h"
#include "langevin.h"
#include "nlistGPU.h"
#include "random.h"
#include "codata.h"
#include "kineticGPU.h"
#include "runtimeKernel.h"
#include "gpu_allocator.hpp"
//void kinetic_terms(SYSTEM*sys, int flag);
//void eval_energyInfo(SYSTEM *sys);

#define PBC(x) x 
#define _DEBUG1(x)  
//alloc memory for constraint data on gpu

void allocCons(SYSTEM * sys, NGLFCONSTRAINT_PARMS *parms)
{
    CHARMMPOT_PARMS *potParms = parms->charmmpot_parms;
    CHARMM_PARMS *charmmParms = potParms->charmmParms;

    //allocate host and device pointers to device memory
    potParms->constraintparms_h = (CONSTRAINTGPU_PARMS *) malloc(sizeof (CONSTRAINTGPU_PARMS));
    potParms->gpu_constraintparms_h = (CONSTRAINTGPU_PARMS*) malloc(sizeof (CONSTRAINTGPU_PARMS));
    gpu_allocator(potParms->gpu_constraintparms, 1);

    //determine necessary size of serialized buffers
    int numPair = 0;
    int numAtom = 0;
    int constraintListSize = 0;
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        RESI_CONN * resiConn = charmmParms->resiConnList[i];
        constraintListSize += resiConn->consListSize;
        for (int j = 0; j < resiConn->consListSize; j++)
        {
            CONSTRAINT *constraint = resiConn->consList[j];
            numPair += constraint->numPair;
            numAtom += constraint->numAtom;
        }
    }

    //alloc host serialized buffers
    CONSTRAINTGPU_PARMS * cParms_h = potParms->constraintparms_h;
    cParms_h->resiConstraintCount = (int *) malloc(sizeof (int)*charmmParms->resiConnSize);
    cParms_h->resiConstraintCountPrefixSum = (int *) malloc(sizeof (int)*charmmParms->resiConnSize);
    cParms_h->constraintPairsPrefixSum = (int *) malloc(sizeof (int)*constraintListSize);
    cParms_h->constraintPairs = (int *) malloc(sizeof (int)*constraintListSize);
    cParms_h->constraintPairsI = (int *) malloc(sizeof (int)*numPair);
    cParms_h->constraintPairsJ = (int *) malloc(sizeof (int)*numPair);
    cParms_h->constraintPairsIindex = (int *) malloc(sizeof (int)*numPair);
    cParms_h->constraintPairsJindex = (int *) malloc(sizeof (int)*numPair);
    cParms_h->constraintPairsDistance = (double *) malloc(sizeof (double)*numPair);
    cParms_h->constraintAtoms = (int *) malloc(sizeof (int)*constraintListSize);
    cParms_h->constraintAtomsPrefixSum = (int *) malloc(sizeof (int)*constraintListSize);
    cParms_h->constraintAtomIDs = (int *) malloc(sizeof (int)*numAtom);

    //alloc gpu buffers
    CONSTRAINTGPU_PARMS *gcp = potParms->gpu_constraintparms_h;
    gpu_allocator(gcp->resiConstraintCount, charmmParms->resiConnSize);
    gpu_allocator(gcp->resiConstraintCountPrefixSum, charmmParms->resiConnSize);
    gpu_allocator(gcp->constraintPairs, constraintListSize);
    gpu_allocator(gcp->constraintPairsPrefixSum, constraintListSize);
    gpu_allocator(gcp->constraintPairsI, numPair);
    gpu_allocator(gcp->constraintPairsJ, numPair);
    gpu_allocator(gcp->constraintPairsIindex, numPair);
    gpu_allocator(gcp->constraintPairsJindex, numPair);
    gpu_allocator(gcp->constraintPairsDistance, numPair);
    gpu_allocator(gcp->constraintAtoms, constraintListSize);
    gpu_allocator(gcp->constraintAtomsPrefixSum, constraintListSize);
    gpu_allocator(gcp->constraintAtomIDs, numAtom);

}

#define PRINT_CDATA(x) 
//serialize all the constraint data 

void serializeCons(SYSTEM * sys, NGLFCONSTRAINT_PARMS *parms)
{

    CHARMMPOT_PARMS *charmmpotParms = parms->charmmpot_parms;
    CHARMM_PARMS *charmmParms = parms->charmmpot_parms->charmmParms;
    CONSTRAINTGPU_PARMS * cParms_h = charmmpotParms->constraintparms_h;

    /*
    //count number of constraints and number of pairs 
    //so we can malloc serial arrays for them
    for (int i = 0; i <charmmParms->resiConnSize; i++)
    {
       RESI_CONN * resiConn = charmmParms->resiConnList[i];
       cParms_h->resiConstraintCountPrefixSum[i] = constraintListSize;
       cParms_h->resiConstraintCount[i] =  resiConn->consListSize;
       constraintListSize += resiConn->consListSize;
       for(int j =0; j< resiConn->consListSize; j++)
       {
          CONSTRAINT *constraint = resiConn->consList[j];
          numPair+= constraint->numPair;
       }
    } 
     */

    //serialize nested data structure
    int ci = 0; //constraint index
    int cpi = 0; //constraint pair index
    int cai = 0; //atomd ID list
    int numAtom2 = 0;
    int numPairs2 = 0; //another counter for determining prefix sum array
    int constraintListSize = 0; //also a counter for making a prefix sum array
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        RESI_CONN * resiConn = charmmParms->resiConnList[i];
        cParms_h->resiConstraintCountPrefixSum[i] = constraintListSize;
        cParms_h->resiConstraintCount[i] = resiConn->consListSize;
        constraintListSize += resiConn->consListSize;
        CONSTRAINT ** consList = resiConn->consList;
        PRINT_CDATA(printf("resi %i   num_constraint %i\n", i, resiConn->consListSize);)
            //iterate over constraint groups in a residue
        for (int j = 0; j < resiConn->consListSize; j++)
        {
            CONSTRAINT *constraint = consList[j];

            cParms_h->constraintPairs[ci] = constraint->numPair;
            cParms_h->constraintPairsPrefixSum[ci] = numPairs2;
            PRINT_CDATA(printf("\tconstraint %i  pairs: %i\n", j, constraint->numPair);)
            numPairs2 += constraint->numPair;
            //iterate over pairs in a constraint group
            for (int k = 0; k < constraint->numPair; k++)
            {
                CONS_PAIR *cpair = constraint->conspairList[k];
                cParms_h->constraintPairsI[cpi] = cpair->atomI;
                cParms_h->constraintPairsJ[cpi] = cpair->atomJ;
                cParms_h->constraintPairsIindex[cpi] = cpair->atomIindex;
                cParms_h->constraintPairsJindex[cpi] = cpair->atomJindex;
                cParms_h->constraintPairsDistance[cpi] = cpair->distance;
                PRINT_CDATA(printf("\t\tpair %i %i \n", cpair->atomI, cpair->atomJ);)
                cpi++;
            }

            cParms_h->constraintAtoms[ci] = constraint->numAtom;
            cParms_h->constraintAtomsPrefixSum[ci] = numAtom2;
            PRINT_CDATA(printf("\tconstraint %i  atoms: %i\n", j, constraint->numAtom);)
            numAtom2 += constraint->numAtom;
            //iterate over atoms in a constraint group
            for (int l = 0; l < constraint->numAtom; l++)
            {
                int atomID = constraint->atomIDList[l];
                cParms_h->constraintAtomIDs[cai] = atomID;
                PRINT_CDATA(printf("\t\t atomID %i\n", cParms_h->constraintAtomIDs[cai]);)
                cai++;
            }
            ci++;
        }
    }
    cParms_h->constraintListSize = constraintListSize;
    cParms_h->numPair = numPairs2;
    cParms_h->numAtom = numAtom2;
    PRINT_CDATA(printf("c list size %i\n", constraintListSize);)

}
#undef PRINT_CDATA

//memcopy serialized cons to gpu

void migrateCons(SYSTEM * sys, NGLFCONSTRAINT_PARMS *parms)
{
    CHARMMPOT_PARMS *potParms = parms->charmmpot_parms;
    CHARMM_PARMS *charmmParms = potParms->charmmParms;
    CONSTRAINTGPU_PARMS * cph = potParms->constraintparms_h;
    CONSTRAINTGPU_PARMS * cpg = potParms->gpu_constraintparms_h;

    int numPair = cph->numPair;
    int numAtom = cph->numAtom;
    int constraintListSize = cph->constraintListSize;
    printf("constraint list size %i pairs %i atoms %i\n", constraintListSize, numPair, numAtom);
    int resiConnSize = charmmParms->resiConnSize;
    gpu_memcpy_host2device(cpg->resiConstraintCount, cph->resiConstraintCount, resiConnSize);
    gpu_memcpy_host2device(cpg->resiConstraintCountPrefixSum, cph->resiConstraintCountPrefixSum, resiConnSize);
    gpu_memcpy_host2device(cpg->constraintPairs, cph->constraintPairs, constraintListSize);
    gpu_memcpy_host2device(cpg->constraintPairsPrefixSum, cph->constraintPairsPrefixSum, constraintListSize);
    gpu_memcpy_host2device(cpg->constraintPairsI, cph->constraintPairsI, numPair);
    gpu_memcpy_host2device(cpg->constraintPairsJ, cph->constraintPairsJ, numPair);
    gpu_memcpy_host2device(cpg->constraintPairsIindex, cph->constraintPairsIindex, numPair);
    gpu_memcpy_host2device(cpg->constraintPairsJindex, cph->constraintPairsJindex, numPair);
    gpu_memcpy_host2device(cpg->constraintPairsDistance, cph->constraintPairsDistance, numPair);
    gpu_memcpy_host2device(cpg->constraintAtoms, cph->constraintAtoms, constraintListSize);
    gpu_memcpy_host2device(cpg->constraintAtomsPrefixSum, cph->constraintAtomsPrefixSum, constraintListSize);
    gpu_memcpy_host2device(cpg->constraintAtomIDs, cph->constraintAtomIDs, numAtom);

    gpu_memcpy_host2device(potParms->gpu_constraintparms, cpg, 1);


}


static RuntimeKernelManager kernelGen;
//static CUfunction kernel1;

void nglfconstraintGPU_parms(SYSTEM * sys, NGLFCONSTRAINT_PARMS *parms)
{
    printf("nglf constraint parms \n");
    if (parms->isotropic == 1) 
        parms->volumeFunc = (void(*)(COLLECTION *collection, STATE *gpu_state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt, int nlocal))changeVolumeGPUisotropic;
     else 
        parms->volumeFunc = (void(*)(COLLECTION *collection, STATE *gpu_state, BOX_STRUCT *box, THREE_SMATRIX *pTensor, double beta, double tau, double dt, int nlocal))changeVolumeGPU;
    charmmResidues(sys, parms->charmmpot_parms);
    allocCons(sys, parms);
    serializeCons(sys, parms);
    migrateCons(sys, parms);
    SETLIST *residueSet = &parms->charmmpot_parms->residueSet;
    int numResidues = residueSet->listSize;
    printf("num Residues %i\n", numResidues);
    printf("nion %i\n", sys->nion);
    parms->starts = (int*) malloc(sizeof (int)*numResidues);
    parms->resID = (int*) malloc(sizeof (int)*sys->nion);
    parms->groupID = (int*) malloc(sizeof (int)*sys->nion);
    parms->consToResStarts = (int*) malloc(sizeof (int)*sys->nion);
    gpu_allocator(parms->startsGPU, numResidues);
    gpu_allocator(parms->resIDGPU, sys->nion);
    gpu_allocator(parms->groupIDGPU, sys->nion);
    gpu_allocator(parms->consToResStartsGPU, sys->nion);

    //runtime kernel
    /*
    std::string fileName = "/g/g19/sundram1/ddcmd/src/constraintKernelBase.cu";
    std::string kernelName = "constraintKernelCustom";
    runtimeKernelSegments segments = kernelGen.parseKernelSegmentFile("/g/g19/sundram1/ddcmd/src/constraintKernelPieces.cu");
    std::string generatedFile = kernelGen.generateKernel(segments, (void*) parms->charmmpot_parms);
    kernelName = "constraintKernelCustomGen";
    kernelGen.compile("outputKernel.cu", kernelName);
    kernelGen.loadToRuntime(kernelName);
    kernelGen.getKernel(&kernel1, "constraintKernelCustomGen");
     */
}

//just for reference
/*
void resMoveConsCopy(int consListSize, CONSTRAINT **consList, double dt, STATE* state, int resRangeStart, int location)
{
   double ( *func) (double dt, double dist2, THREE_VECTOR rab,THREE_VECTOR vab);
   //if (location == FRONT_TIMESTEP) func = frontFunc; 
   //if (location == BACK_TIMESTEP ) func = backFunc; 
   const double tol=1.0e-8; //unit less
   const int maxit=500;

   for(int i=0; i<consListSize; i++)
   {
      CONSTRAINT* constraint = consList[i]; 
      //printf("group=%d numPair=%d\n",i,constraint->numPair); 
      int numAtm = constraint->numAtom; // total number of atom
      double   rMass[numAtm]; 
      THREE_VECTOR r[numAtm],v[numAtm]; 
      for(int j = 0; j < numAtm; ++j)
      {
         int index=resRangeStart+constraint->atomIDList[j];

         rMass[j]=1.0/((ATOMTYPE_PARMS *) (state->species[index]->parm))->mass;
         VSET(r[j],state->rx[index],state->ry[index],state->rz[index]);
         VSET(v[j],state->vx[index],state->vy[index],state->vz[index]);
      }
      int numPair = constraint->numPair; 
      THREE_VECTOR rab[numPair];
      for (int ab = 0; ab < numPair; ++ab) 
      {
         CONS_PAIR* consPair=constraint->conspairList[ab];
         int a=consPair->atomIindex;
         int b=consPair->atomJindex;
         VOP2(rab[ab],=,r[a],-,r[b]); 
         nearestImage(&rab[ab].x, &rab[ab].y, &rab[ab].z);
      }
      //int it=0;
      double errMax = 0.0;    
      for (int it=0;it<maxit;it++) //start the iterative loop
      {
         errMax = 0.0;    
         for (int ab = 0; ab < numPair; ++ab) 
         {
            CONS_PAIR* consPair=constraint->conspairList[ab];
            int a=consPair->atomIindex;
            int b=consPair->atomJindex;
            double dist2=SQ(consPair->distance);
            THREE_VECTOR vab; 
            VOP2(vab,=,v[a],-,v[b]); 
            double rma=rMass[a];
            double rmb=rMass[b];                       

            double rvab = func(dt,dist2, rab[ab], vab)/dist2 ; 

            double gab=-rvab/(rma+rmb);   //units: mass/time)
            double err = fabs(rvab*dt);
            if (err > errMax ) errMax = err;  
            VSVOP(v[a],+=,(rma*gab),*,rab[ab]); 
            VSVOP(v[b],-=,(rmb*gab),*,rab[ab]); 
            
           //if (err < tol) printf("%d: group=%d ab= %d numPair=%d it=%d %e %e\n",getRank(0),i,ab,numPair,it,err, gab); 
         } 
         if (errMax < tol ) 
         {
            break;   
         }
      } //end of iterative loop
    //  if(it == maxit) printf("%d resMove: too many contraint iterations.\n",getRank(0));

      for(int j = 0; j < numAtm; ++j) // Store the new values
      {
         //printf("vfx %f vy %f vz %f\n", v[j].x, v[j].y, v[j].z);
         int index=resRangeStart+constraint->atomIDList[j];
         state->vx[index]=v[j].x;
         state->vy[index]=v[j].y;
         state->vz[index]=v[j].z;   
      }            
   }
}
 */
// we can make these host device functions TODO

__device__ double frontFunc(double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab)
{
    THREE_VECTOR pab = rab;
    VSVOP(pab, +=, dt, *, vab);
    double pab2 = VSQ(pab);
    double rvab = (pab2 - dist2) / (2 * dt);
    //             p2-d2 = (r + v*dt)^2 - d2  = 0.5*(r2-d2)/dt + r.v +0.5*dt*v2 

    return rvab;
}

__device__ double backFunc(double dt, double dist2, THREE_VECTOR rab, THREE_VECTOR vab)
{
    double rvab = DOT(rab, vab);
    return rvab;
}

template <typename Tin, typename Vec>
__device__ Tin frontFuncOpt(Tin dt, Tin dt2inv, Tin dist2, Vec rab, Vec vab)
{
    Vec pab = rab;
    VSVOP(pab, +=, dt, *, vab);
    Tin pab2 = VSQ(pab);
    Tin rvab = (pab2 - dist2) * dt2inv;
    //             p2-d2 = (r + v*dt)^2 - d2  = 0.5*(r2-d2)/dt + r.v +0.5*dt*v2 

    return rvab;
}



//same as above, just to have same argument structure as frontFuncOpt

template <typename Tin, typename Vec>
__device__ Tin backFuncOpt(Tin dt, Tin dt2inv, Tin dist2, Vec rab, Vec vab)
{
    Tin rvab = DOT(rab, vab);
    return rvab;
}




//we need a count of the #residues in system
//and then a count of the number of constraint in the system
//the number of warps launched = #constraints
//we need a resID 

//strategy 1
//each warp takes care of a constraint  
//16 threads per warp, CHOL and PAP have only 1/1 constraint, each w/ 6/7 pairs
//and assign 1 thread/pair
//FOR KRAS, assuming there are >32 pairs, the warp may have to do a 2nd or 3rd  iteration of work
//so each thread will get incremented by warpsize

//strategy 2
//each warp owns a residue
//CHOL and PAP are the same as strategy 1
//for KRAS the warp will have to iterator over pairs, and possibly over constraints as well
//this seems expensive, but it makes scheduling harder
// because then we need a buffer that can map a warp id -> constraint id

// so let's go with strategy 1
//this means we need an array of size(residues in system)
//that contain each residue's resID
//int resID = resIDMap[wid]

//then we need to get the corresponding number of constrains for that residue
//int numConstraints = resiConstraintCount[resID]

//backInBoxGPUO3(THREE_MATRIX *h1, THREE_MATRIX *h1i, double *x, double *y, double *z)
#if (__CUDA_ARCHH__ == 70)
#define SYNCH_GROUP() __syncwarp()
#define BALLOT(x,y) __ballot_sync(x,y)
#else
#define SYNCH_GROUP() __syncthreads()
#define BALLOT(x,y) __ballot(y); __syncthreads()
#endif
#define WARP_SIZE 32

__global__ void constraintKernelOld(STATE *state, THREE_MATRIX *h1, CharmmBONDGPU_PARMS *bparms, CONSTRAINTGPU_PARMS *parms,
                                    gid_type *label, int *species, int * resIDmap, int *residueStarts,
                                    int *consToResStarts, int *groupIDs, double *mass, double dt, int isFrontStep)
{
    //int pid = blockIdx.x*blockDim.x+threadIdx.x;
    int wid = blockIdx.x; //warp id. assuming 32 threads per warp
    int tid = threadIdx.x;
    //asumme numAtm<WARP_SIZE

    __shared__ double rMass[WARP_SIZE]; //numAtm
    __shared__ THREE_VECTOR r[WARP_SIZE]; //numAtm
    __shared__ THREE_VECTOR v[WARP_SIZE]; //numAtm
    __shared__ THREE_VECTOR rab[WARP_SIZE]; //numPair
    __shared__ double gamma[WARP_SIZE];
    const double tol = 1.0e-8; //unit less
    const int maxit = 500;

    //grab pointers 
    //int *numAtomsInResidue = bparms->numAtomsInResidue;

    //every  particle needs a resID
    int resID = wid;
    //int groupId = wid;
    //every resID needs a start index
    int start = residueStarts[resID];
    //int start = consToResStarts[wid];
    //printf("resID %i index %i  tid %i \n", resID, index, tid);
    //choose which residue type to operate on

    int resTypeID = resIDmap[start];
    //if (tid>=numAtomsInResidue[resTypeID]) return;

    int *resiConstraintCount = parms->resiConstraintCount;
    int *resiConstraintCountPrefixSum = parms->resiConstraintCountPrefixSum;
    //int *numPairsPrefixSum = parms->resiConstraintCountPrefixSum;
    int *constraintPairs = parms->constraintPairs;
    int *constraintPairsPrefixSum = parms->constraintPairsPrefixSum;
    //int *constraintPairsI = parms->constraintPairsI;
    //int *constraintPairsJ = parms->constraintPairsJ;
    int *constraintPairsIindex = parms->constraintPairsIindex;
    int *constraintPairsJindex = parms->constraintPairsJindex;
    double *constraintPairsDistance = parms->constraintPairsDistance;
    int *constraintAtoms = parms->constraintAtoms;
    int *constraintAtomsPrefixSum = parms->constraintAtomsPrefixSum;
    int *constraintAtomIDs = parms->constraintAtomIDs;



    //int groupID=groupIDs[wid];
    //int resTypeID = resID[index]; 
    //choose which residue to operate on
    //int numAtm = numAtomsInResidue[resTypeID]; // total number of atom  
    int resi_constraints_start = resiConstraintCountPrefixSum[resTypeID];
    int resi_constraints_end = resiConstraintCountPrefixSum[resTypeID] + resiConstraintCount[resTypeID];

    int cid = 0;
    //iterate over residue's constraints
    for (int constraintID = resi_constraints_start; constraintID < resi_constraints_end; constraintID++)
    {
        //int constraintID = resi_constraints_start+groupID;
        //if(tid==0 && wid==0)printf("wid %i constraintID %i start %i group %i\n", wid, constraintID, resi_constraints_start, groupID);
        //int constraintID = resi_constraints_start + (wid - resiConstraintCountPrefixSum[resTypeID]); 
        //double rMass[];
        //THREE_VECTOR r[numAtm],v[numAtm];  
        int atomIDs_start = constraintAtomsPrefixSum[constraintID];
        int pairIDs_start = constraintPairsPrefixSum[constraintID];
        int numAtm = constraintAtoms[constraintID];
        int numPair = constraintPairs[constraintID];
        if (numPair == 0)
        {
            continue;
        }//continue;}
        cid++;
        //load constraint group information and positions/velocities of that constraint group into shared memory
        if (tid < numAtm)
        {
            int atomID = constraintAtomIDs[atomIDs_start + tid];
            int index = start + atomID;
            double rxi = state->rx[index];
            ;
            double ryi = state->ry[index];
            double rzi = state->rz[index];
            PBC(preduceGPUO3(h1, &rxi, &ryi, &rzi);)
            //iterate over atms in current residue
            r[tid].x = rxi;
            r[tid].y = ryi;
            r[tid].z = rzi;
            v[tid].x = state->vx[index];
            v[tid].y = state->vy[index];
            v[tid].z = state->vz[index];
            int speciesID = species[index];
            rMass[tid] = mass[speciesID];
        }
        SYNCH_GROUP();


        //for KRAS we will have an additonal loop here that iterates of groups of 32 pairs
        //use this for now  
        //for (int p =0; p<numPair; p++)
        int ab = tid;
        if (ab < numPair)
        {
            int pairID = constraintPairsPrefixSum[constraintID] + ab;
            pairID = pairIDs_start + ab;
            int a = constraintPairsIindex[pairID];
            int b = constraintPairsJindex[pairID];
            //double distance = constraintPairsDistance[pairID];
            //double dist2 = SQ(distance);
            VOP2(rab[ab], =, r[a], -, r[b]);
            PBC(preduceGPUO3(h1, &rab[ab].x, &rab[ab].y, &rab[ab].z);)
                //if(resID==0)printf("start resID %i type %i  ab %i a %i b %i va %f vb %f rab %f\n", resID, resTypeID, ab, a, b,v[a].x, v[b].x, rab[ab].x );
        }//return;
        double errMax = 0.0;
        int it = 0;
        //for (int ab = 0; ab < numPair; ++ab) gamma[ab]=0.0;
        if (ab < numPair)
        {
            gamma[ab] = 0.0;
        }
        SYNCH_GROUP();
        //return;
        if (tid == 0)
        {
            for (; it < maxit; it++) //start the iterative loop
            {
                errMax = 0.0;
                for (int ab = 0; ab < numPair; ab++)
                {
                    int pairID = constraintPairsPrefixSum[constraintID] + ab;
                    int a = constraintPairsIindex[pairID];
                    int b = constraintPairsJindex[pairID];
                    double distance = constraintPairsDistance[pairID];
                    double dist2 = SQ(distance);
                    THREE_VECTOR vab;
                    VOP2(vab, =, v[a], -, v[b]);
                    double rma = rMass[a];
                    double rmb = rMass[b];
                    //this shouldn't cause warp divergence because all threads
                    //in kernel will call the same function
                    //this is thus a good candidate for templating TODO
                    double rvab = 0;
                    if (isFrontStep)
                    {
                        rvab = frontFunc(dt, dist2, rab[ab], vab) / dist2;
                    }
                    else
                    {
                        rvab = backFunc(dt, dist2, rab[ab], vab) / dist2;
                    }
                    double gab = -rvab / (rma + rmb); //units: mass/time)
                    double err = fabs(rvab * dt);
                    //if (err > errMax ) errMax = err;  
                    errMax = fmax(err, errMax);
                    VSVOP(v[a], +=, (rma * gab), *, rab[ab]);
                    VSVOP(v[b], -=, (rmb * gab), *, rab[ab]);
                    gamma[ab] += gab;
                } //numPair
                if (tid == 0 && errMax < tol)
                {
                    break;
                }

                //SYNCH_GROUP();        
                /*
                //32 threads/warp in hex
                #define WARP_MASK 0xffff
                int ballot =1; 
                if (tid<numAtm){
                  ballot = (errMax < tol);
                }
                unsigned ballot_result = BALLOT(WARP_MASK, ballot);
                if (!(ballot_result==0xffff))
                {
                   break;
                }
                #undef WARP_MASK
                 */
            } //maxit
        }
        // if (tid==0 && resID ==1115) printf("wid %i tid %i resID %i cid %i range %i %i  resType %i numPair %i iter %i\n",wid, tid, resID, cid, resi_constraints_start, resi_constraints_end, resTypeID,numPair, it);
        SYNCH_GROUP();

        //for(int j = 0; j < numAtm; ++j) // Store the new values
        if (tid < numAtm)
        {
            //int index=resRangeStart+constraint->atomIDList[j];
            int atomID = constraintAtomIDs[atomIDs_start + tid];
            int index = start + atomID;

            state->vx[index] = v[tid].x;
            state->vy[index] = v[tid].y;
            state->vz[index] = v[tid].z;
        }
        atomIDs_start += numAtm;
        pairIDs_start += numPair;
        SYNCH_GROUP();
    } //constraints in resieu

}


// signature for all valid template params
typedef double(*constraint_function)(double dt, double dt2inv, double dist2, THREE_VECTOR rab, THREE_VECTOR vab);
//typedef float(*constraint_function32)(float dt,float dt2inv, float dist2, THREE_VECTORF rab,THREE_VECTORF vab); 

/*
kernel contains 3 optimizations over constraintKernelOld
   1)use of shared memory to store intermediate values in iterative solver
   2)use of template function to call frontFunc or backFunc, instead of an if statement
   3)ability to use FP32 to store intermediate values in iterative solver instead of FT64
     this is controlled using FP, Vec, and constraint_function template parameters
 */
template <constraint_function cons_func, typename FP, typename Vec>
__global__ void constraintKernel(STATE *state, THREE_MATRIX *h1, CharmmBONDGPU_PARMS *bparms, CONSTRAINTGPU_PARMS *parms,
                                 gid_type *label, int *species, int * resIDmap, int *residueStarts,
                                 int *consToResStarts, int *groupIDs, double *mass, double dtt, int nConstraint, int isFrontStep)
{

    //int pid = blockIdx.x*blockDim.x+threadIdx.x;
    int wid = blockIdx.x; //warp id. assuming 32 threads per warp
    int tid = threadIdx.x;
    if (wid >= nConstraint)return;
    //asumme numAtm<WARP_SIZE
    FP dt = dtt;
    __shared__ FP rMass[WARP_SIZE]; //numAtm
    __shared__ Vec r[WARP_SIZE]; //numAtm
    __shared__ Vec v[WARP_SIZE]; //numAtm
    __shared__ Vec rab[WARP_SIZE]; //numPair
    __shared__ FP rMassPairInv[WARP_SIZE]; //numPair
    __shared__ FP gamma[WARP_SIZE];
    //__shared__ int pairIDs[WARP_SIZE];
    __shared__ int pairsIindexs[WARP_SIZE];
    __shared__ int pairsJindexs[WARP_SIZE];
    __shared__ FP distances2[WARP_SIZE];
    __shared__ FP distances2inv[WARP_SIZE];
    const FP tol = 1.0e-12; //unit less
    const int maxit = 500;

    //grab pointers 
    //int *numAtomsInResidue = bparms->numAtomsInResidue;

    //every  particle needs a resID
    int resID = wid;

    //every resID needs a start index
    int start = residueStarts[resID];
    start = consToResStarts[wid];
    //choose which residue type to operate on

    int resTypeID = resIDmap[start];
    //if (tid>=numAtomsInResidue[resTypeID]) return;

    //int *resiConstraintCount = parms->resiConstraintCount;
    int *resiConstraintCountPrefixSum = parms->resiConstraintCountPrefixSum;
    //int *numPairsPrefixSum = parms->resiConstraintCountPrefixSum;
    int *constraintPairs = parms->constraintPairs;
    int *constraintPairsPrefixSum = parms->constraintPairsPrefixSum;
    //int *constraintPairsI = parms->constraintPairsI;
    //int *constraintPairsJ = parms->constraintPairsJ;
    int *constraintPairsIindex = parms->constraintPairsIindex;
    int *constraintPairsJindex = parms->constraintPairsJindex;
    double *constraintPairsDistance = parms->constraintPairsDistance;
    int *constraintAtoms = parms->constraintAtoms;
    int *constraintAtomsPrefixSum = parms->constraintAtomsPrefixSum;
    int *constraintAtomIDs = parms->constraintAtomIDs;
    int groupID = groupIDs[wid];

    //choose which residue to operate on
    int resi_constraints_start = resiConstraintCountPrefixSum[resTypeID];
    //int resi_constraints_end = resiConstraintCountPrefixSum[resTypeID]+resiConstraintCount[resTypeID];

    int cid = 0;
    //iterate over residue's constraints
    //for (int constraintID=resi_constraints_start; constraintID<resi_constraints_end; constraintID++)
    {
        int constraintID = resi_constraints_start + groupID;
        //if(tid==0 && wid==1004)printf("wid %i constraintID %i start %i group %i\n", wid, constraintID, resi_constraints_start, groupID);
        int atomIDs_start = constraintAtomsPrefixSum[constraintID];
        int pairIDs_start = constraintPairsPrefixSum[constraintID];
        int numAtm = constraintAtoms[constraintID];
        int numPair = constraintPairs[constraintID];
        if (numPair == 0)
        {
            return;
        }//continue;}
        cid++;
        if (tid < numAtm)
        {
            int atomID = constraintAtomIDs[atomIDs_start + tid];
            int index = start + atomID;
            FP rxi = state->rx[index];
            ;
            FP ryi = state->ry[index];
            FP rzi = state->rz[index];
            FP h1xx = h1->xx;
            FP h1yy = h1->yy;
            FP h1zz = h1->zz;
            PBC(preduceGPUO3x(h1xx, h1yy, h1zz, rxi, ryi, rzi);)
            //iterate over atms in current residue
            r[tid].x = rxi;
            r[tid].y = ryi;
            r[tid].z = rzi;
            v[tid].x = state->vx[index];
            v[tid].y = state->vy[index];
            v[tid].z = state->vz[index];
            int speciesID = species[index];
            rMass[tid] = mass[speciesID];
        }
        SYNCH_GROUP();


        //for KRAS we will have an additonal loop here that iterates of groups of 32 pairs
        //use this for now  
        //for (int p =0; p<numPair; p++)
        FP dt2inv = 1.0 / (2 * dt);
        int ab = tid;
        if (ab < numPair)
        {
            int pairID = constraintPairsPrefixSum[constraintID] + ab;
            pairID = pairIDs_start + ab;
            int a = constraintPairsIindex[pairID];
            int b = constraintPairsJindex[pairID];
            FP distance = constraintPairsDistance[pairID];
            FP dist2 = SQ(distance);
            VOP2(rab[ab], =, r[a], -, r[b]);
            PBC(preduceGPUO3(h1, &rab[ab].x, &rab[ab].y, &rab[ab].z);)
            //if(resID==0)printf("start resID %i type %i  ab %i a %i b %i va %f vb %f rab %f\n", resID, resTypeID, ab, a, b,v[a].x, v[b].x, rab[ab].x );
            //pairIDs[ab]=pairID;
            pairsIindexs[ab] = a;
            pairsJindexs[ab] = b;
            distances2[ab] = dist2;
            distances2inv[ab] = 1.0 / dist2;
            rMassPairInv[ab] = 1.0 / (rMass[a] + rMass[b]);

        }//return;
        SYNCH_GROUP();
        FP errMax = 0.0;
        int it = 0;
        if (ab < numPair)
        {
            gamma[ab] = 0.0;
        }
        SYNCH_GROUP();
        if (tid == 0)
        {
            for (; it < maxit; it++) //start the iterative loop
            {
                errMax = 0.0;
                for (int ab = 0; ab < numPair; ab++)
                {
                    //int pairID = constraintPairsPrefixSum[constraintID]+ab;
                    //int a = constraintPairsIindex[pairID];
                    //int b = constraintPairsJindex[pairID];
                    //FP distance = constraintPairsDistance[pairID];
                    //int pairID = pairIDs[ab];
                    int a = pairsIindexs[ab];
                    int b = pairsJindexs[ab];
                    FP dist2 = distances2[ab];
                    //FP dist2 = SQ(distance);
                    FP dist2inv = distances2inv[ab];
                    Vec vab;
                    VOP2(vab, =, v[a], -, v[b]);
                    //_DEBUG1(if(resID==0 ) printf("ab %i a %i b %i vax %f vay %f vaz %f\n", ab, a, b, vab.x, vab.y, vab.z);)
                    FP rma = rMass[a];
                    FP rmb = rMass[b];
                    FP rma_rmb_inv = rMassPairInv[ab];
                    //this shouldn't cause warp divergence because all threads
                    //in kernel will call the same function
                    //this is thus a good candidate for templating TODO
                    FP rvab = 0;
                    rvab = cons_func(dt, dt2inv, dist2, rab[ab], vab) * dist2inv;
                    //FP gab=-rvab/(rma+rmb);   //units: mass/time)
                    FP gab = -rvab*rma_rmb_inv; //units: mass/time)
                    FP err = fabs(rvab * dt);
                    errMax = fmax(err, errMax);
                    VSVOP(v[a], +=, (rma * gab), *, rab[ab]);
                    VSVOP(v[b], -=, (rmb * gab), *, rab[ab]);
                    gamma[ab] += gab;
                } //numPair
                if (tid == 0 && errMax < tol)
                {
                    break;
                }

            } //maxit
        }
        SYNCH_GROUP();

        if (tid < numAtm)
        {
            int atomID = constraintAtomIDs[atomIDs_start + tid];
            int index = start + atomID;

            state->vx[index] = v[tid].x;
            state->vy[index] = v[tid].y;
            state->vz[index] = v[tid].z;
        }
        atomIDs_start += numAtm;
        pairIDs_start += numPair;
        SYNCH_GROUP();
    } //constraints in resieu

}

__global__ void constraintKernelOld(STATE *state, THREE_MATRIX *h1, CharmmBONDGPU_PARMS *bparms, CONSTRAINTGPU_PARMS *parms,
                                    gid_type *label, int *species, int * resIDmap, int *residueStarts,
                                    int *consToResStarts, int *groupIDs, double *mass, double dt, int nConstraint, int isFrontStep)
{
    //int pid = blockIdx.x*blockDim.x+threadIdx.x;
    int wid = blockIdx.x; //warp id. assuming 32 threads per warp
    int tid = threadIdx.x;
    if (wid >= nConstraint)return;
    //asumme numAtm<WARP_SIZE

    __shared__ double rMass[WARP_SIZE]; //numAtm
    __shared__ THREE_VECTOR r[WARP_SIZE]; //numAtm
    __shared__ THREE_VECTOR v[WARP_SIZE]; //numAtm
    __shared__ THREE_VECTOR rab[WARP_SIZE]; //numPair
    __shared__ double gamma[WARP_SIZE];
    const double tol = 1.0e-12; //unit less
    const int maxit = 500;

    //grab pointers 
    //int *numAtomsInResidue = bparms->numAtomsInResidue;

    //every  particle needs a resID
    int resID = wid;

    //every resID needs a start index
    int start = residueStarts[resID];
    start = consToResStarts[wid];
    //choose which residue type to operate on

    int resTypeID = resIDmap[start];
    //if (tid>=numAtomsInResidue[resTypeID]) return;

    //int *resiConstraintCount = parms->resiConstraintCount;
    int *resiConstraintCountPrefixSum = parms->resiConstraintCountPrefixSum;
    //int *numPairsPrefixSum = parms->resiConstraintCountPrefixSum;
    int *constraintPairs = parms->constraintPairs;
    int *constraintPairsPrefixSum = parms->constraintPairsPrefixSum;
    //int *constraintPairsI = parms->constraintPairsI;
    //int *constraintPairsJ = parms->constraintPairsJ;
    int *constraintPairsIindex = parms->constraintPairsIindex;
    int *constraintPairsJindex = parms->constraintPairsJindex;
    double *constraintPairsDistance = parms->constraintPairsDistance;
    int *constraintAtoms = parms->constraintAtoms;
    int *constraintAtomsPrefixSum = parms->constraintAtomsPrefixSum;
    int *constraintAtomIDs = parms->constraintAtomIDs;
    int groupID = groupIDs[wid];

    //choose which residue to operate on
    int resi_constraints_start = resiConstraintCountPrefixSum[resTypeID];
    //int resi_constraints_end = resiConstraintCountPrefixSum[resTypeID]+resiConstraintCount[resTypeID];

    int cid = 0;
    //iterate over residue's constraints
    //for (int constraintID=resi_constraints_start; constraintID<resi_constraints_end; constraintID++)
    {
        int constraintID = resi_constraints_start + groupID;
        //if(tid==0 && wid==1004)printf("wid %i constraintID %i start %i group %i\n", wid, constraintID, resi_constraints_start, groupID);
        int atomIDs_start = constraintAtomsPrefixSum[constraintID];
        int pairIDs_start = constraintPairsPrefixSum[constraintID];
        int numAtm = constraintAtoms[constraintID];
        int numPair = constraintPairs[constraintID];
        if (numPair == 0)
        {
            return;
        }//continue;}
        cid++;
        _DEBUG1(
                if (tid == 0 && resID == 1115)
                printf("tid %i wid %i start %i constraintID %i end %i  count %i atm %i pair %i\n", tid, wid, resi_constraints_start, constraintID, resi_constraints_end, resiConstraintCount[wid], numAtm, numPair);
                )

            if (tid < numAtm)
            {
                int atomID = constraintAtomIDs[atomIDs_start + tid];
                int index = start + atomID;
                double rxi = state->rx[index];
                ;
                double ryi = state->ry[index];
                double rzi = state->rz[index];
                PBC(preduceGPUO3(h1, &rxi, &ryi, &rzi);)
                //iterate over atms in current residue
                r[tid].x = rxi;
                r[tid].y = ryi;
                r[tid].z = rzi;
                v[tid].x = state->vx[index];
                v[tid].y = state->vy[index];
                v[tid].z = state->vz[index];
                _DEBUG1(if (resID == 1115)printf("\tresID %i resTypeID %i tid %i index %i atomID %i ais %i numatm %i vx %f\n", resID, resTypeID, tid, index, atomID, atomIDs_start, numAtm, v[tid].x);)
                    int speciesID = species[index];
                rMass[tid] = mass[speciesID];
            }
        SYNCH_GROUP();


        //for KRAS we will have an additonal loop here that iterates of groups of 32 pairs
        //use this for now  
        //for (int p =0; p<numPair; p++)
        int ab = tid;
        if (ab < numPair)
        {
            int pairID = constraintPairsPrefixSum[constraintID] + ab;
            pairID = pairIDs_start + ab;
            int a = constraintPairsIindex[pairID];
            int b = constraintPairsJindex[pairID];
            //double distance = constraintPairsDistance[pairID];
            //double dist2 = SQ(distance);
            VOP2(rab[ab], =, r[a], -, r[b]);
            PBC(preduceGPUO3(h1, &rab[ab].x, &rab[ab].y, &rab[ab].z);)
                //if(resID==0)printf("start resID %i type %i  ab %i a %i b %i va %f vb %f rab %f\n", resID, resTypeID, ab, a, b,v[a].x, v[b].x, rab[ab].x );
        }//return;
        double errMax = 0.0;
        int it = 0;
        if (ab < numPair)
        {
            gamma[ab] = 0.0;
        }
        SYNCH_GROUP();
        if (tid == 0)
        {
            for (; it < maxit; it++) //start the iterative loop
            {
                errMax = 0.0;
                for (int ab = 0; ab < numPair; ab++)
                {
                    int pairID = constraintPairsPrefixSum[constraintID] + ab;
                    int a = constraintPairsIindex[pairID];
                    int b = constraintPairsJindex[pairID];
                    double distance = constraintPairsDistance[pairID];
                    double dist2 = SQ(distance);
                    THREE_VECTOR vab;
                    VOP2(vab, =, v[a], -, v[b]);
                    //_DEBUG1(if(resID==0 ) printf("ab %i a %i b %i vax %f vay %f vaz %f\n", ab, a, b, vab.x, vab.y, vab.z);)
                    double rma = rMass[a];
                    double rmb = rMass[b];
                    //this shouldn't cause warp divergence because all threads
                    //in kernel will call the same function
                    //this is thus a good candidate for templating TODO
                    double rvab = 0;
                    if (isFrontStep)
                    {
                        rvab = frontFunc(dt, dist2, rab[ab], vab) / dist2;
                    }
                    else
                    {
                        rvab = backFunc(dt, dist2, rab[ab], vab) / dist2;
                    }
                    _DEBUG1(if (resID == 1115)printf("\tit %i a %i b %i va %f vb %f rvab %f\n", it, a, b, v[a].x, v[b].x, rvab);)
                        double gab = -rvab / (rma + rmb); //units: mass/time)
                    double err = fabs(rvab * dt);
                    errMax = fmax(err, errMax);
                    VSVOP(v[a], +=, (rma * gab), *, rab[ab]);
                    VSVOP(v[b], -=, (rmb * gab), *, rab[ab]);
                    gamma[ab] += gab;
                } //numPair
                if (tid == 0 && errMax < tol)
                {
                    break;
                }

            } //maxit
        }
        // if (tid==0 && resID ==1115) printf("wid %i tid %i resID %i cid %i range %i %i  resType %i numPair %i iter %i\n",wid, tid, resID, cid, resi_constraints_start, resi_constraints_end, resTypeID,numPair, it);
        SYNCH_GROUP();

        if (tid < numAtm)
        {
            int atomID = constraintAtomIDs[atomIDs_start + tid];
            int index = start + atomID;

            state->vx[index] = v[tid].x;
            state->vy[index] = v[tid].y;
            state->vz[index] = v[tid].z;
        }
        atomIDs_start += numAtm;
        pairIDs_start += numPair;
        SYNCH_GROUP();
    } //constraints in resieu

}

void updateStateAliasesGPU(SYSTEM *sys, unsigned *nlocal, double **rx, double **ry, double **rz, double **vx, double **vy, double **vz, double **fx, double **fy, double **fz, SPECIES ***species, gid_type **label)
{
    STATE *state = sys->collection->state;
    *nlocal = sys->nlocal;
    *rx = state->rx;
    *ry = state->ry;
    *rz = state->rz; // The SYSTEM and STATE might change during the call to ddcenergy
    *vx = state->vx;
    *vy = state->vy;
    *vz = state->vz; // (i.e., we might reassign particles to new tasks) so we need to
    *fx = state->fx;
    *fy = state->fy;
    *fz = state->fz; // update all of the aliases we use.
    *label = state->label;
    *species = state->species;
}

/*
void getResiStarts(int *starts, int *resID, STATE *state, int numResidues, RESRANGE *resRange2)
{
   int constraint=0;
   for(int i =0; i<numResidues;i++)
   {
     starts[i]=resRange2[i].start;
   }
}
 */

int generateSendResStartMap(SYSTEM *sys, NGLFCONSTRAINT_PARMS*p)
{
    STATE * state = sys->collection->state;
    //calculate particle->residue mapping arrays for kernle
    //this only needs to be done once per redomain
    SETLIST *residueSet = &p->charmmpot_parms->residueSet;
    LISTNODE* residueList = residueSet->list;
    unsigned start = 0;
    int numConstraints = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        p->starts[i] = start;
        start += residueList[i].size;
        RESI_CONN *resiConn = residueList[i].resiConn;
        for (int j = 0; j < resiConn->consListSize; j++)
        {
            if (resiConn->consList[j]->numPair > 0)
            {
                numConstraints++;
            }
        }
    }

    //calculate a map from each constraint group to it's corresponding
    //residue start
    int consI = 0;
    for (int i = 0; i < residueSet->listSize; i++)
    {
        RESI_CONN *resiConn = residueList[i].resiConn;
        for (int j = 0; j < resiConn->consListSize; j++)
        {
            if (resiConn->consList[j]->numPair > 0)
            {
                //printf("consI %i name %s numPair %i\n", consI, resiConn->resName, resiConn->consList[j]->numPair);
                p->consToResStarts[consI] = p->starts[i];
                p->groupID[consI] = j;
                consI++;
            }
        }
    }
    //printf("num constraints %i\n", numConstraints);
    int numResidues = residueSet->listSize;
    //getResiStarts(p->starts,p->resID,  state, numResidues, resRange2);
    gpu_memcpy_host2device(p->startsGPU, p->starts, numResidues);
    gpu_memcpy_host2device(p->resIDGPU, p->charmmpot_parms->bondparms_h->resID, state->nlocal);
    gpu_memcpy_host2device(p->consToResStartsGPU, p->consToResStarts, numConstraints);
    gpu_memcpy_host2device(p->groupIDGPU, p->groupID, numConstraints);
    return consI;
}

__global__ void backInBoxKernel(THREE_MATRIX *h1, THREE_MATRIX *h1i, STATE *state, int nlocal)
{

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid < nlocal)
    {
        backInBoxGPUO3(h1, h1i, &state->rx[pid], &state->ry[pid], &state->rz[pid]);
    }
}

void launchRuntimeConstraintKernel(COLLECTION *collection, STATE *gsh, NGLFCONSTRAINT_PARMS *p, CHARMMPOT_PARMS *cp_parms, int numResidues,
        const CUfunction& kernelFunc, double dt) 
{
    GPUNLIST *gnlist = collection->gnlist;

    CUdeviceptr d_rx = (CUdeviceptr) gsh->rx;
    CUdeviceptr d_ry = (CUdeviceptr) gsh->ry;
    CUdeviceptr d_rz = (CUdeviceptr) gsh->rz;

    CUdeviceptr d_vx = (CUdeviceptr) gsh->vx;
    CUdeviceptr d_vy = (CUdeviceptr) gsh->vy;
    CUdeviceptr d_vz = (CUdeviceptr) gsh->vz;

    CUdeviceptr d_h1 = (CUdeviceptr) gnlist->hmat_g;
    CUdeviceptr d_cpparms = (CUdeviceptr) cp_parms->gpu_constraintparms;

    CUdeviceptr d_species = (CUdeviceptr) gnlist->species_g;
    CUdeviceptr d_resIDmap = (CUdeviceptr) p->resIDGPU;
    CUdeviceptr d_resStarts = (CUdeviceptr) p->startsGPU;

    CUdeviceptr d_consToResStarts = (CUdeviceptr) p->consToResStartsGPU;
    CUdeviceptr d_groupID = (CUdeviceptr) p->groupIDGPU;
    CUdeviceptr d_mass = (CUdeviceptr) gnlist->mass_g;
    //CUdeviceptr d_dt = (CUdeviceptr) dt;
    CUdeviceptr d_nConstraint = (CUdeviceptr) numResidues;
    CUdeviceptr dtg;
    cuMemAlloc(&dtg, sizeof (double));
    cuMemcpyHtoD(dtg, &dt, sizeof (double));

    int blockSizeCons = WARP_SIZE;
    int gridSizeCons = numResidues;
    //printf("dt 1 %f \n", dt);
    void *args[] = {&d_rx, &d_ry, &d_rz,
        &d_vx, &d_vy, &d_vz,
        &d_h1, &d_cpparms,
        &d_species, &d_resIDmap, &d_resStarts,
        &d_consToResStarts, &d_groupID, &d_mass, &dtg, &d_nConstraint};



    cuLaunchKernel(kernelFunc,
                   gridSizeCons, 1, 1, // grid dim
                   blockSizeCons, 1, 1, // block dim
                   0, NULL, // shared mem and stream
                   args, // arguments
                   0);
    cuCtxSynchronize();
}




static int numResidues;

void nglfconstraintGPU(DDC*ddc, SIMULATE*simulate, NGLFCONSTRAINT_PARMS*p) 
{
    PUSH_RANGE("timestep", 4);
    double dt = simulate->dt;
    double time = simulate->time;
    SYSTEM* sys = simulate->system;
    GPUNLIST *gnlist = sys->collection->gnlist;
    STATE* state = sys->collection->state;
    STATE* gsh = sys->collection->gpustate_h;
    COLLECTION *collection = sys->collection;
    unsigned nlocal = state->nlocal;
    double *rx, *ry, *rz, *vx, *vy, *vz, *fx, *fy, *fz;
    SPECIES** species;
    LONG64* label;

    calcMolecularPressuresGPU(sys);

    calcTensorWithGPUVirials(sys, &sys->energyInfo);

    THREE_SMATRIX pTensor = (sys->energyInfo.virial);

    int N = sys->moleculeClass->nMoleculesGlobal;
    double vol = box_get_volume(NULL);
    double T = p->T;
    pTensor.xx += N * kB*T;
    pTensor.yy += N * kB*T;
    pTensor.zz += N * kB*T;
    SMATNORM(pTensor, vol);
    pTensor.xx -= p->P0;
    pTensor.yy -= p->P0;
    pTensor.zz -= p->P0;

    //printf("gpu pTensor %.15f %.15f %.15f\n\n", pTensor.xx, pTensor.yy, pTensor.zz);
    p->volumeFunc(collection, collection->gpustate, sys->box, &pTensor, p->beta, p->tauBarostat, dt, sys->nlocal);
    //printf("ho %f %f %f \n", sys->box->h0.xx, sys->box->h0.yy, sys->box->h0.zz);



    int blockSize = min(64, state->nlocal);
    int gridSize = ceil((float) state->nlocal / blockSize);
    //freeVelocityUpdateGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->mass_g, collection->gpustate, state->nlocal, .5*dt);
    freeVelocityUpdateInterleavedGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->listIdsg, gnlist->results, gnlist->mass_g, collection->gpustate, state->nlocal, .5 * dt);

    PUSH_RANGE("regen map", 2);
    if (ddc->lastUpdate == sys->loop)
    { // re-assign map upon domain changes
        numResidues = generateSendResStartMap(sys, p);
        THREE_MATRIX *hinvPtr = &(sys->box->hinv);
        gpu_memcpy_host2device(gnlist->hmati_g, hinvPtr, 1);
        backInBoxKernel<<<gridSize, blockSize>>>(gnlist->hmat_g, gnlist->hmati_g, collection->gpustate, state->nlocal);
        CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }
    POP_RANGE();

    CHARMMPOT_PARMS *cp_parms = p->charmmpot_parms;
    int blockSizeCons = WARP_SIZE;
    int gridSizeCons = numResidues;
    if (gridSizeCons > 0) {
        //constraintKernelOld<<<gridSizeCons,blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
        constraintKernel<frontFuncOpt, double, THREE_VECTOR><<<gridSizeCons, blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
                cp_parms->gpu_constraintparms, gsh->label, gnlist->species_g, p->resIDGPU,
                p->startsGPU, p->consToResStartsGPU, p->groupIDGPU, gnlist->mass_g, dt, numResidues, 1);
        //launchRuntimeConstraintKernel(collection, gsh, p, cp_parms, numResidues, kernel1, dt);
        CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }

    updateStateAliasesGPU(sys, &nlocal, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    freePositionUpdateGPU<<<gridSize, blockSize>>>(collection->gpustate, state->nlocal, dt);
    CUDA_SAFE_CALL(cudaPeekAtLastError();)
    backInBoxKernel<<<gridSize, blockSize>>>(gnlist->hmat_g, gnlist->hmati_g, collection->gpustate, state->nlocal);
    //double time_plus_dt = time + dt; 


    ddc->update = 0;
    time += dt; // positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt 
    simulate->time = sys->time = time;
    simulate->loop++;
    sys->loop = simulate->loop;
    //for (int kk=0;kk<sys->ngroup;kk++) sys->group[kk]->Update1(sys->group[kk],-1,state,time,0.5*dt);
    if (ddcenergy(ddc, sys, 0) != 0) return;


    updateStateAliasesGPU(sys, &nlocal, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    //freeVelocityUpdateGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->mass_g, collection->gpustate,nlocal, .5*dt);
    freeVelocityUpdateInterleavedGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->listIdsg, gnlist->results, gnlist->mass_g, collection->gpustate, state->nlocal, .5 * dt);
    CUDA_SAFE_CALL(cudaPeekAtLastError();)
            //constraintKernelOld<<<gridSizeCons,blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
    if (gridSizeCons > 0) {
        constraintKernel<backFuncOpt, double, THREE_VECTOR><<<gridSizeCons, blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
                cp_parms->gpu_constraintparms, gsh->label, gnlist->species_g, p->resIDGPU,
                p->startsGPU, p->consToResStartsGPU, p->groupIDGPU, gnlist->mass_g, dt, numResidues, 0);
        CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }
    //kineticGPU(sys, 1);


    /*errorCheck(ddc->domain_id, simulate->loop, state, sys->energyInfo, p, datafile); */
    simulate->time = sys->time = time;
    POP_RANGE();
}

void nglfconstraintGPULangevin(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS*p)
{
    PUSH_RANGE("timestep", 4);
    double dt = simulate->dt;
    double time = simulate->time;
    SYSTEM* sys = simulate->system;
    GPUNLIST *gnlist = sys->collection->gnlist;
    STATE* state = sys->collection->state;
    STATE* gsh = sys->collection->gpustate_h;
    COLLECTION *collection = sys->collection;
    unsigned nlocal = state->nlocal;
    double *rx, *ry, *rz, *vx, *vy, *vz, *fx, *fy, *fz;
    SPECIES** species;
    LONG64* label;

    calcMolecularPressuresGPU(sys);

    calcTensorWithGPUVirials(sys, &sys->energyInfo);

    THREE_SMATRIX pTensor = (sys->energyInfo.virial);

    //printf("pTensor %f %f %f\n", pTensor.xx, pTensor.yy, pTensor.zz);
    int N = sys->moleculeClass->nMoleculesGlobal;
    double vol = box_get_volume(NULL);
    double T = p->T;
    pTensor.xx += N * kB*T;
    pTensor.yy += N * kB*T;
    pTensor.zz += N * kB*T;
    SMATNORM(pTensor, vol);
    pTensor.xx -= p->P0;
    pTensor.yy -= p->P0;
    pTensor.zz -= p->P0;

    //printf("gpu pTensor %.15f %.15f %.15f\n\n", pTensor.xx, pTensor.yy, pTensor.zz);

    p->volumeFunc(collection, collection->gpustate, sys->box, &pTensor, p->beta, p->tauBarostat, dt, sys->nlocal);

    //get random numbers
    PUSH_RANGE("random gen", 3);
    unsigned nrandom = nlocal;
    if (nrandom % 2 != 0) nrandom++;

    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomx, nrandom, 0.0, 1.0);
    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomy, nrandom, 0.0, 1.0);
    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomz, nrandom, 0.0, 1.0);
    POP_RANGE();

    GROUP* group0 = state->group[0];
    LANGEVIN_PARMS *lp = (LANGEVIN_PARMS*) group0->parm;
    THREE_VECTOR vcm = lp->vcm;
    double tau = lp->tau[0];
    double kBT = lp->kBT[0];

    int blockSize = min(64, state->nlocal);
    //int gridSize = (state->nlocal+blockSize-1)/blockSize;
    int gridSize = ceil((float) state->nlocal / blockSize);
    //freeVelocityUpdateGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->mass_g, collection->gpustate, state->nlocal, .5*dt);
    langevin_velocityUpdateGPUInterleaved_FrontTimestep<<<gridSize, blockSize>>>(gnlist->species_g, collection->gpustate, gnlist->results, tau, kBT, gnlist->mass_g, state->nlocal, .5 * dt, vcm.x, vcm.y, vcm.z, gnlist->gpuRandomx, gnlist->gpuRandomy, gnlist->gpuRandomz);
    if (ddc->lastUpdate == sys->loop) // re-assign map upon domain changes
        numResidues = generateSendResStartMap(sys, p);

    if (sys->loop % 20 == 0)
    { // re-assign map upon domain changes
        THREE_MATRIX *hinvPtr = &(sys->box->hinv);
        gpu_memcpy_host2device(gnlist->hmati_g, hinvPtr, 1);
        backInBoxKernel<<<gridSize, blockSize>>>(gnlist->hmat_g, gnlist->hmati_g, collection->gpustate, state->nlocal);
        CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }

    CHARMMPOT_PARMS *cp_parms = p->charmmpot_parms;
    int blockSizeCons = WARP_SIZE;
    int gridSizeCons = numResidues;
    if (gridSizeCons > 0) {
        //constraintKernelOld<<<gridSizeCons,blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
        constraintKernel<frontFuncOpt, double, THREE_VECTOR><<<gridSizeCons, blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
                cp_parms->gpu_constraintparms, gsh->label, gnlist->species_g, p->resIDGPU,
                p->startsGPU, p->consToResStartsGPU, p->groupIDGPU, gnlist->mass_g, dt, numResidues, 1);
        CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }
    updateStateAliasesGPU(sys, &nlocal, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    freePositionUpdateGPU<<<gridSize, blockSize>>>(collection->gpustate, state->nlocal, dt);

    CUDA_SAFE_CALL(cudaPeekAtLastError();)
    //double time_plus_dt = time + dt; 

    ddc->update = 0;
    time += dt; // positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt 
    simulate->time = sys->time = time;
    simulate->loop++;
    sys->loop = simulate->loop;

    if (ddcenergy(ddc, sys, 0) != 0) return;

    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomx, nrandom, 0.0, 1.0);
    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomy, nrandom, 0.0, 1.0);
    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomz, nrandom, 0.0, 1.0);

    vcm = lp->vcm;
    updateStateAliasesGPU(sys, &nlocal, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    langevin_velocityUpdateGPUInterleaved_BackTimestep<<<gridSize, blockSize>>>(gnlist->species_g, collection->gpustate, gnlist->results, tau, kBT, gnlist->mass_g, state->nlocal, .5 * dt, vcm.x, vcm.y, vcm.z, gnlist->gpuRandomx, gnlist->gpuRandomy, gnlist->gpuRandomz);
    //freeVelocityUpdateGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->mass_g, collection->gpustate,nlocal, .5*dt);
    CUDA_SAFE_CALL(cudaPeekAtLastError();)
    if (gridSizeCons > 0) {
        //constraintKernelOld<<<gridSizeCons,blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
        constraintKernel<backFuncOpt, double, THREE_VECTOR><<<gridSizeCons, blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
                cp_parms->gpu_constraintparms, gsh->label, gnlist->species_g, p->resIDGPU,
                p->startsGPU, p->consToResStartsGPU, p->groupIDGPU, gnlist->mass_g, dt, numResidues, 0);
        CUDA_SAFE_CALL(cudaPeekAtLastError();)
    }
    //kineticGPU(sys, 1);
    simulate->time = sys->time = time;
    POP_RANGE();

}

/*
FOR DEBUGING ONLY
Same as nglfconstraintGPULangevin, but generates the random numbers on CPU
and then memcopies them to GPU. Useful for checking energies against CPU code
 */
void nglfconstraintGPULangevinLCG64(DDC *ddc, SIMULATE *simulate, NGLFCONSTRAINT_PARMS*p) 
{
    PUSH_RANGE("timestep", 4);
    double dt = simulate->dt;
    double time = simulate->time;
    SYSTEM* sys = simulate->system;
    GPUNLIST *gnlist = sys->collection->gnlist;
    STATE* state = sys->collection->state;
    STATE* gsh = sys->collection->gpustate_h;

    COLLECTION *collection = sys->collection;
    unsigned nlocal = state->nlocal;
    double *rx, *ry, *rz, *vx, *vy, *vz, *fx, *fy, *fz;
    SPECIES** species;
    LONG64* label;

    calcMolecularPressuresGPU(sys);

    calcTensorWithGPUVirials(sys, &sys->energyInfo);

    THREE_SMATRIX pTensor = (sys->energyInfo.virial);

    //printf("pTensor %f %f %f\n", pTensor.xx, pTensor.yy, pTensor.zz);
    int N = sys->moleculeClass->nMoleculesGlobal;
    double vol = box_get_volume(NULL);
    double T = p->T;
    pTensor.xx += N * kB*T;
    pTensor.yy += N * kB*T;
    pTensor.zz += N * kB*T;
    SMATNORM(pTensor, vol);
    pTensor.xx -= p->P0;
    pTensor.yy -= p->P0;
    pTensor.zz -= p->P0;

    //printf("gpu pTensor %.15f %.15f %.15f\n\n", pTensor.xx, pTensor.yy, pTensor.zz);

    p->volumeFunc(collection, collection->gpustate, sys->box, &pTensor, p->beta, p->tauBarostat, dt, sys->nlocal);
    //printf("ho %f %f %f \n", sys->box->h0.xx, sys->box->h0.yy, sys->box->h0.zz);

    double randNumsx[sys->nlocal];
    double randNumsy[sys->nlocal];
    double randNumsz[sys->nlocal];

    //double *randNumsxg, *randNumsyg, *randNumszg;
    /*
    gpu_allocator(randNumsxg,  sys->nlocal);
    gpu_allocator(randNumsyg,  sys->nlocal);
    gpu_allocator(randNumszg,  sys->nlocal);
     */
    THREE_VECTOR g;
    //generate random numbers

    for (int i = 0; i < sys->nlocal; i++)
    {

        GROUP* group = state->group[0];
        LANGEVIN_PARMS *p = (LANGEVIN_PARMS*) group->parm;
        RANDOM *random = p->random;
        void *randomParms = random_getParms(random, i);
        g = gasdev3d(random, randomParms);
        randNumsx[i] = g.x;
        randNumsy[i] = g.y;
        randNumsz[i] = g.z;
        //if (i<10) printf("rand %i %f %f %f\n", i, g.x, g.y, g.z);
    }

    /*  
    gpu_memcpy_host2device(randNumsxg,  randNumsx,  sys->nlocal);
    gpu_memcpy_host2device(randNumsyg,  randNumsy,  sys->nlocal);
    gpu_memcpy_host2device(randNumszg,  randNumsz,  sys->nlocal);
     */
    double *randNumsxPtr = randNumsx;
    double *randNumsyPtr = randNumsy;
    double *randNumszPtr = randNumsz;
    gpu_memcpy_host2device(gnlist->gpuRandomx, randNumsxPtr, sys->nlocal);
    gpu_memcpy_host2device(gnlist->gpuRandomy, randNumsyPtr, sys->nlocal);
    gpu_memcpy_host2device(gnlist->gpuRandomz, randNumszPtr, sys->nlocal);


    GROUP* group0 = state->group[0];
    LANGEVIN_PARMS *lp = (LANGEVIN_PARMS*) group0->parm;
    THREE_VECTOR vcm = lp->vcm;
    double tau = lp->tau[0];
    double kBT = lp->kBT[0];

    int blockSize = min(64, state->nlocal);
    int gridSize = ceil((float) state->nlocal / blockSize);
    //freeVelocityUpdateGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->mass_g, collection->gpustate, state->nlocal, .5*dt);
    langevin_velocityUpdateGPUInterleaved_FrontTimestep<<<gridSize, blockSize>>>(gnlist->species_g, collection->gpustate, gnlist->results, tau, kBT, gnlist->mass_g, state->nlocal, .5 * dt, vcm.x, vcm.y, vcm.z, gnlist->gpuRandomx, gnlist->gpuRandomy, gnlist->gpuRandomz);
    if (ddc->lastUpdate == sys->loop) { // re-assign map upon domain changes
        numResidues = generateSendResStartMap(sys, p);
    }
    CHARMMPOT_PARMS *cp_parms = p->charmmpot_parms;

    int blockSizeCons = WARP_SIZE;
    int gridSizeCons = numResidues;

    constraintKernelOld<<<gridSizeCons, blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
            cp_parms->gpu_constraintparms, gsh->label, gnlist->species_g, p->resIDGPU,
            p->startsGPU, p->consToResStartsGPU, p->groupIDGPU, gnlist->mass_g, dt, numResidues, 1);

    //launchRuntimeConstraintKernel(collection, gsh, p, cp_parms, kernel1, dt);

    CUDA_SAFE_CALL(cudaPeekAtLastError();)
    updateStateAliasesGPU(sys, &nlocal, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    freePositionUpdateGPU<<<gridSize, blockSize>>>(collection->gpustate, state->nlocal, dt);

    //generate random numbers
    for (int i = 0; i < sys->nlocal; i++)
    {

        GROUP* group = state->group[0];
        LANGEVIN_PARMS *p = (LANGEVIN_PARMS*) group->parm;
        RANDOM *random = p->random;
        void *randomParms = random_getParms(random, i);
        g = gasdev3d(random, randomParms);
        randNumsx[i] = g.x;
        randNumsy[i] = g.y;
        randNumsz[i] = g.z;
        //if (i<10) printf("rand %i %f %f %f\n", i, g.x, g.y, g.z);

    }

    /*
    gpu_memcpy_host2device(randNumsxg,  randNumsx,  sys->nlocal);
    gpu_memcpy_host2device(randNumsyg,  randNumsy,  sys->nlocal);
    gpu_memcpy_host2device(randNumszg,  randNumsz,  sys->nlocal);
     */

    gpu_memcpy_host2device(gnlist->gpuRandomx, randNumsxPtr, sys->nlocal);
    gpu_memcpy_host2device(gnlist->gpuRandomy, randNumsyPtr, sys->nlocal);
    gpu_memcpy_host2device(gnlist->gpuRandomz, randNumszPtr, sys->nlocal);


    CUDA_SAFE_CALL(cudaPeekAtLastError();)

    //double time_plus_dt = time + dt; 

    ddc->update = 0;
    time += dt; // positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt 
    simulate->time = sys->time = time;
    simulate->loop++;
    sys->loop = simulate->loop;

    if (ddcenergy(ddc, sys, 0) != 0) return;
    /*
    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomx, nrandom, 0.0, 1.0);
    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomy, nrandom, 0.0, 1.0);
    curandGenerateNormalDouble(gnlist->gpu_types->gen, gnlist->gpuRandomz, nrandom, 0.0, 1.0);
     */
    //printf("flag\n");
    //cudaDeviceSynchronize();
    vcm = lp->vcm;
    updateStateAliasesGPU(sys, &nlocal, &rx, &ry, &rz, &vx, &vy, &vz, &fx, &fy, &fz, &species, &label);
    langevin_velocityUpdateGPUInterleaved_BackTimestep<<<gridSize, blockSize>>>(gnlist->species_g, collection->gpustate, gnlist->results, tau, kBT, gnlist->mass_g, state->nlocal, .5 * dt, vcm.x, vcm.y, vcm.z, gnlist->gpuRandomx, gnlist->gpuRandomy, gnlist->gpuRandomz);
    //freeVelocityUpdateGPU<<<gridSize, blockSize>>>(gnlist->species_g, gnlist->mass_g, collection->gpustate,nlocal, .5*dt);
    CUDA_SAFE_CALL(cudaPeekAtLastError();)
    constraintKernelOld<<<gridSizeCons, blockSizeCons>>>(collection->gpustate, gnlist->hmat_g, cp_parms->gpu_bondparms,
            cp_parms->gpu_constraintparms, gsh->label, gnlist->species_g, p->resIDGPU,
            p->startsGPU, p->consToResStartsGPU, p->groupIDGPU, gnlist->mass_g, dt, numResidues, 0);
    CUDA_SAFE_CALL(cudaPeekAtLastError();)


    simulate->time = sys->time = time;
    POP_RANGE();

    //cudaFree(randNumsxg);
    //cudaFree(randNumsyg);
    //cudaFree(randNumszg);


}

/*
void nglfconstraintGPULangevinOld(DDC*ddc, SIMULATE*simulate, NGLFCONSTRAINT_PARMS*p)
{
        double dt = simulate->dt;
        double time = simulate->time;
        SYSTEM* sys = simulate->system;
   COLLECTION *collection = sys->collection;
        STATE* state = sys->collection->state;
        STATE* gsh = sys->collection->gpustate_h;
   unsigned nlocal=state->nlocal; 
        double *rx,*ry,*rz,*vx,*vy,*vz,*fx,*fy,*fz;
        SPECIES** species; 
        LONG64* label ; 


   int blockSize = min(64,state->nlocal) ;
   int gridSize = ceil((float) state->nlocal/blockSize); 
   unsigned nrandom = nlocal;
   if(nrandom % 2 != 0) nrandom++;
   

   curandGenerateNormalDouble(collection->gpu_types->gen, collection->gpuRandomx, nrandom, 0.0, 1.0 );

   curandGenerateNormalDouble(collection->gpu_types->gen, collection->gpuRandomy, nrandom, 0.0, 1.0);
   curandGenerateNormalDouble(collection->gpu_types->gen, collection->gpuRandomz, nrandom, 0.0, 1.0);

   GROUP *group = state->group[0];
   LANGEVIN_PARMS *lp=(LANGEVIN_PARMS *)group->parm;

   THREE_VECTOR vcm = lp->vcm;
   double tau = lp->tau[0];
   double kBT = lp->kBT[0];
   //Front Timestep
  // langevin_velocityUpdateGPU_FrontTimestep<<<gridSize, blockSize>>>(collection->species_g, collection->gpustate, tau, kBT, collection->mass_g, nlocal, 0.5*dt, vcm.x, vcm.y, vcm.z, collection->gpuRandomx, collection->gpuRandomy, collection->gpuRandomz);


   //freeVelocityUpdateGPU<<<gridSize, blockSize>>>(collection->species_g, collection->mass_g, collection->gpustate, state->nlocal, .5*dt);
 
  
   //calculate particle->residue mapping arrays for kernle
   //this only needs to be done once
   
   //charmmResidues(sys, p->charmmpot_parms);
   SETLIST *residueSet=&p->charmmpot_parms->residueSet;
   LISTNODE* residueList=residueSet->list;
   RESRANGE resRange2[residueSet->listSize];  // Need to put on heap
   unsigned start =0;
   for (int i =0 ; i <residueSet->listSize; i ++)
   {
      resRange2[i].start = start;
      start+=residueList[i].size;
      resRange2[i].end = start;		
   }
   int numResidues = residueSet->listSize;
   int *starts = (int*) malloc(sizeof(int)*numResidues);
   int *resID = (int*) malloc(sizeof(int)*state->nion);
 
   //getResiStarts(starts,resID, state, numResidues, resRange2);
   //printf("num residues %i\n", numResidues);
   int *startsGPU;
   int *resIDGPU;
   gpu_allocator(startsGPU, numResidues);
   gpu_allocator(resIDGPU, state->nion);
   gpu_memcpy_host2device(startsGPU,  starts, numResidues);
   gpu_memcpy_host2device(resIDGPU,  p->charmmpot_parms->bondparms_h->resID, state->nion);
   //CUDA_SAFE_CALL(cudaMemcpy(collection->hmat_g, &(sys->box->h0), sizeof(THREE_MATRIX), cudaMemcpyHostToDevice);)
   CHARMMPOT_PARMS *cp_parms = p->charmmpot_parms;
 
   //printf("first constraint update\n");
   int blockSizeCons = 32;
   int gridSizeCons = residueSet->listSize; 
   constraintKernel<<<gridSizeCons,blockSizeCons>>>(collection->gpustate, collection->hmat_g, cp_parms->gpu_bondparms,
                               cp_parms->gpu_constraintparms, gsh->label, collection->species_g, resIDGPU, startsGPU, 
                               p->consToResStartsGPU,  p->groupIDGPU,collection->mass_g, dt, 1);
   CUDA_SAFE_CALL(cudaPeekAtLastError();)
   updateStateAliasesGPU(sys,&nlocal,&rx,&ry,&rz,&vx,&vy,&vz,&fx,&fy,&fz,&species,&label);
 

   freePositionUpdateGPU<<<gridSize, blockSize>>>(collection->gpustate,state->nlocal, dt);
   sendHostState(sys, nlocal);
   double time_plus_dt = time + dt; 

   ddc->update = 0;
   time += dt;                                     // positions, box (volume, h0,hinv) , and forces at  t = n*dt + dt 
   simulate->time=sys->time=time; 
   simulate->loop++;
   sys->loop = simulate->loop;
   if (ddcenergy(ddc, sys, 0) != 0) return;
   updateStateAliasesGPU(sys,&nlocal,&rx,&ry,&rz,&vx,&vy,&vz,&fx,&fy,&fz,&species,&label);
   //sendForceVelocityToGPU(sys, nlocal);
//   freeVelocityUpdateGPU<<<gridSize, blockSize>>>(collection->species_g, collection->mass_g, collection->gpustate,nlocal, .5*dt);
   //get new random numbers and do back timestep
      curandGenerateNormalDouble(collection->gpu_types->gen, collection->gpuRandomx, nrandom, 0.0, 1.0);
   curandGenerateNormalDouble(collection->gpu_types->gen, collection->gpuRandomy, nrandom, 0.0, 1.0);
   curandGenerateNormalDouble(collection->gpu_types->gen, collection->gpuRandomz, nrandom, 0.0, 1.0);
//   langevin_velocityUpdateGPU_BackTimestep<<<gridSize, blockSize>>>(collection->species_g, collection->gpustate, tau, kBT, collection->mass_g, nlocal, 0.5*dt, vcm.x, vcm.y, vcm.z, collection->gpuRandomx, collection->gpuRandomy, collection->gpuRandomz);

   constraintKernel<<<gridSizeCons,blockSizeCons>>>(collection->gpustate, collection->hmat_g, cp_parms->gpu_bondparms,
                               cp_parms->gpu_constraintparms, gsh->label, collection->species_g, resIDGPU, 
                              startsGPU, p->consToResStartsGPU, p->groupIDGPU, collection->mass_g, dt, 0);
   CUDA_SAFE_CALL(cudaPeekAtLastError();)

   simulate->time = sys->time = time;
   free(starts); free(resID);
}
 */
/* Local Variables: */
/* tab-width: 3 */
/* End: */
