//#include <thrust/sort.h>
//#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <nvToolsExt.h>
#include "ddcMalloc.h"
#include "pairProcessGPU.h"
#include "nlistGPU.h"
#include "cudaUtils.h"
#include "gpuMemUtils.h"
#include "prefixScan.h"
#include "pair.h"
#include "energyInfo.h"
#include "minMaxGPU.h"
#include <stdio.h>
#include <stdlib.h>
#include "molecularPressureGPU.h"
#if (USE_GPU ==1)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#endif
#include <algorithm>
#include <sys/time.h>
#include <cmath>
#include <float.h>
#include "gpu_allocator.hpp"
#include "accelerator.h"
//pair process
#define SHARED 128
#define SHARED_EVAL 32
#define SHARED_EVALp 32
#define TPP_EVAL 1
#define TPP_EVALp 1
#define  CHECK_1(x)     //first layer of Debugging, summary of binning data
#define CHECK_2(x)    //second layer of debugging, print all binning data
#define CHECK_3(x)   //not used rn
#define CHECK(x)   //final layer of debugging, print all particle data
#define CHECKB (x) 
#define CHECK2(x) 
#define GPROF(x)  //used when doing profiling w/o nvprof
#define TPP //threads per particle

/*
This file is responsible for doing particle binning operations
on the GPU. Most of the binning operations happen in a C function
binParticles() which calls a series of GPU kernels 

1. calculate the destination bin of every particle 
   (where each bin edge is of length rcut) ; "assignBinIDs()" Kernel
2. calculates the number of particles that will reside in every bin 
   and the 'offset' each particle relative to the start of the bin in memory
   "getBinCounts" kernel
   The kernels in 3. and 4. live in the binSort() funcion, which calls a series of kernels
3. does prefix sum on the particlecount-per-bin array calculated in 2.
4. adds the 'offset' from 2. and the particle-counts perfix-sum value from 3. to 
   get new particle position in bin order
Finally, the particles are "permuted" into binned order by a function called
permuteParticles(), that stores the permute particle data in seperate gpu buffers
ie rxbg, rybg, rzbg for particle positions
 */

//based on the dimensions (lens) of each uniformly sized bin, 
//and the number of said bins in each dimension (nbins)
//this function determines the bin index for each particle
//in each dimension (bx,by,bz). It then flattens these indeces into 1 index

__global__ void assignBinIDs(unsigned *bin,
                             double *rx, double *ry, double *rz, double * lens,
                             int *nbins, double *mins, int n)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= n)
    {
        return;
    }
    //divide the particle position (relative to simulation box origin)
    // by bin length to determine bin index in each dimension
    int bx = nbins[0]*(rx[pid] - mins[0]) / (lens[0]);
    int by = nbins[1]*(ry[pid] - mins[1]) / (lens[1]);
    int bz = nbins[2]*(rz[pid] - mins[2]) / (lens[2]);

    bx = min(max(bx, 0), nbins[0] - 1);
    by = min(max(by, 0), nbins[1] - 1);
    bz = min(max(bz, 0), nbins[2] - 1);
    bin[pid] = (unsigned) (nbins[0]*(nbins[1] * bz + by) + bx); //flatten indeces
}


//initializes r_back. simply stores array index in each array position

__global__ void fillBackPointers(int *r_backg, int n)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= n)
    {
        return;
    }
    r_backg[pid] = pid;
}

//simply uses index mapping from particle to binned particle index
//stored in r_back to copy particles into new 'binned' positions

__global__ void permuteParticlesKernel(int *new_ids, double *rx, double *ry, double *rz, int *r_back,
                                       double *rxbg, double *rybg, double *rzbg, unsigned *bins, unsigned *binsb, int *speciesg,
                                       int *speciesbg, double *charge_g, double *charge_bg, int n)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= n)
    {
        return;
    }
    int new_id = new_ids[pid];
    rxbg[new_id] = rx[pid];
    rybg[new_id] = ry[pid];
    rzbg[new_id] = rz[pid];
    r_back[new_id] = pid;
    binsb[new_id] = bins[pid];
    speciesbg[new_id] = speciesg[pid];
    charge_bg[new_id] = charge_g[pid];
}

typedef struct ewaldgpu
{
    double eps;
    double sig;
} EWALDGPU;

/*
void permuteParticles(int *new_ids, double *rx, double *ry, double *rz, int *r_back,
                              double *rxbg,double *rybg,double *rzbg,unsigned *bins,unsigned *binsb, int n)
{
   
   int blocksize = SHARED;
   int gridsize = ceil((float)n/blocksize);
   permuteParticlesKernel<<<gridsize, blocksize>>>(new_ids, rx, ry, rz, r_back, rxbg, rybg, rzbg, bins, binsb,n);
}
 */
void permuteParticles2(COLLECTION *collection, int n) {

    int blocksize = SHARED;
    int gridsize = ceil((float) n / blocksize);
    GPUNLIST *gnlist=collection->gnlist;
    STATE *gsh = collection->gpustate_h;
    permuteParticlesKernel<<<gridsize, blocksize>>>(gnlist->listIdsg, gsh->rx, gsh->ry, gsh->rz, gnlist->r_backg, gnlist->rxbg, gnlist->rybg, gnlist->rzbg, gnlist->binsg, gnlist->binsbg, gnlist->species_g, gnlist->species_bg, gnlist->charge_g, gnlist->charge_bg, n);

   /*
            int nion = n;
            int *r_back = (int*) malloc(nion * sizeof (int));
            gpu_memcpy_device2host(r_back, gnlist->listIdsg, nion);
    for (int i = 0; i < nion; i++) {
        for (int j = i+1; j < nion; j++) {
       		if(r_back[i]==r_back[j]) printf("%d  %d listid %d  %d\n", i, j, r_back[i], r_back[j]);
        }
    }
    */
    /*
    int bsize=256;
    int gsize=ceil((float) nion / bsize);
    printf("Debug Pair orig\n\n");

    debugPair<<<gsize, bsize>>>(gsh->rx, gsh->ry, gsh->rz, nion);

    printf("Debug Pair new\n\n");
    debugPair<<<gsize, bsize>>>(gnlist->rxbg, gnlist->rybg, gnlist->rzbg, nion);
     */
}

void permuteParticlesCharmm(COLLECTION *collection, CHARMMPOT_PARMS *parms, int n) {
    //
    CharmmLJGPU_PARMS *gpu_parms = parms->gpu_ljparms_h;
    int blocksize = SHARED;
    int gridsize = ceil((float) n / blocksize);
    GPUNLIST *gnlist=collection->gnlist;
    STATE *gsh = collection->gpustate_h;
    permuteParticlesKernel<<<gridsize, blocksize>>>(gnlist->listIdsg, gsh->rx, gsh->ry, gsh->rz, gnlist->r_backg, gnlist->rxbg, gnlist->rybg, gnlist->rzbg, gnlist->binsg, gnlist->binsbg, gpu_parms->cg_species_index, gpu_parms->cg_species_index_b, gnlist->charge_g, gnlist->charge_bg, n);
}
//Lennard Jones potential

__device__ thrust::pair<double, double> Kernel(LJGPU parms, double r)
{
    double sigma_r = parms.sig * (1 / r);
    double eps = parms.eps;
    double valid = r <= parms.rcut;
    double shift = parms.shift;
    double s2 = sigma_r*sigma_r;
    double s4 = s2*s2;
    double s6 = s4*s2;
    double s12 = s6*s6;

    double dvdr_over_r = valid * 4.0 * eps * (-12.0 * s12 + 6.0 * s6) / (r * r);
    double e = valid * (4.0 * eps * (s12 - s6) + shift);
    thrust::pair<double, double> p(e, dvdr_over_r);
    return p;
}

//Lennard Jones pair.c potential without divisions. The thrust stuff needs to go

__device__ thrust::pair<double, double> NKernel2(LJGPU parms, double r, double rinv, double rinv2)
{
    double sigma_r = parms.sig*rinv;
    double eps = parms.eps;
    double valid = r < parms.rcut;
    double shift = parms.shift;
    double s2 = sigma_r*sigma_r;
    double s4 = s2*s2;
    double s6 = s4*s2;
    double s12 = s6*s6;

    double dvdr_over_r = valid * 4.0 * eps * (-12.0 * s12 + 6.0 * s6) * rinv2;
    double e = valid * (4.0 * eps * (s12 - s6) + shift);
    thrust::pair<double, double> p(e, dvdr_over_r);
    return p;
}




//Ewald SR potential

__device__ thrust::pair<double, double> Kernel(EWALDGPU parms, double r)
{
    //double sigma = reinterpret_cast<double> parms[0];
    //double eps = reinterpret_cast<double> parms[1];
    thrust::pair<double, double> p(0, 0);
    return p;
}


//this function uses a non-neighbor list technique to evaluate lennard jones
//it is not optimized
//I also don't know if this works

template <typename K>
__global__ void processPairNaiveOld(int *rback, double *rxbg, double *rybg, double *rzbg,
                                    double *e, double *fx, double *fy, double *fz,
                                    unsigned *bin, int *nbins, int *nbrIdsx, int *nbrIdsy, int *nbrIdsz,
                                    int *binHeads, int* numNbrs, int *binCounts, K * parms, int n)
{
    int check = 4;
    //1 thread spawned per particle atm, this could change
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    //shouldn't have more threads than particles
    if (pid >= n)
    {
        return;
    }

    int ii = rback[pid]; //use back pointer to get original particle index
    K p = *parms;
    double f = 0;
    int binsTot = nbins[0] * nbins[1] * nbins[2];
    int local_bin = bin[pid];
    CHECK(
          if (pid == check)
    {
          printf("local bin %i\n", local_bin);
    }
          if (pid == check)
    {
          printf("bin_tot %i\n", binsTot);
    }
          )
    //get local particle assigned to this thread
    double rxl = rxbg[pid];
    double ryl = rybg[pid];
    double rzl = rzbg[pid];

    int loc_bin = local_bin;
    //if bin is on the boundary, some neighboring bins wont exist. skip those bins
    //this shouldn't result in branch diversion if each warp hits same condition
    int xybins = nbins[0] * nbins[1];
    int loc_bin_z = local_bin / (xybins);

    int loc_bin_xy = (loc_bin - loc_bin_z * (xybins));

    int loc_bin_y = loc_bin_xy / (nbins[0]);
    int loc_bin_x = loc_bin_xy - loc_bin_y * (nbins[0]);

    CHECK(if (pid == check)
    {
          printf("z y x %i %i %i\n", numNbrs[0], numNbrs[1], numNbrs[2]);
    })

    double ee = 0;
    double ffx = 0;
    double ffy = 0;
    double ffz = 0;


    //loop over 27 possible neighboring bins
    for (int i = 0; i < numNbrs[0]; i++)
    {
        for (int j = 0; j < numNbrs[1]; j++)
        {
            for (int k = 0; k < numNbrs[2]; k++)
            {
                //if bin is on the boundary, some neighboring bins wont exist. skip those bins
                //this shouldn't result in branch diversion if each warp hits same condition
                int nbr_bin_z = loc_bin_z + nbrIdsz[k];
                int nbr_bin_y = loc_bin_y + nbrIdsy[j];
                int nbr_bin_x = loc_bin_x + nbrIdsx[i];
                int nbr_bin = nbr_bin_z * nbins[0] * nbins[1] + nbr_bin_y * nbins[0] + nbr_bin_x;
                CHECK(if (pid == check)
                {
                      printf("nbr_bin %i, %i %i %i\n", nbr_bin, nbr_bin_x, nbr_bin_y, nbr_bin_z); })
                if (nbr_bin_x < 0 || nbr_bin_x >= nbins[0] || nbr_bin_y < 0 || nbr_bin_y >= nbins[1] || nbr_bin_z < 0 || nbr_bin_z >= nbins[2])
                {
                    continue;
                }
                CHECK(if (pid == check)
                {
                      printf("g_bin %i, %i %i %i\n", nbr_bin, nbr_bin_x, nbr_bin_y, nbr_bin_z); })
                //now we have to loop over remote particles in neighboring bin.  
                //we process them in "chunks" of size=blocksize     
                //int chunks = ceil((float) binCounts[nbr_bin]/blockDim.x);
                //for (int c = 0;c<chunks; c++)
                //{
                //this is where we would import to shared memory
                //int z = binHeads[nbr_bin]+c*blockDim.x+threadIdx.x
                //sdata_rx[tid] = rxbg[z]; ...etc etc
                CHECK(if (pid == check)
                {
                      printf("binC binH %i %i  \n\n", binCounts[nbr_bin], binHeads[nbr_bin]);
                })
                int bin_end = binHeads[nbr_bin] + binCounts[nbr_bin];
                //binHeads points to start of bin. binCounts is length of bin
                for (int j = binHeads[nbr_bin]; j < binHeads[nbr_bin] + binCounts[nbr_bin]; j++)
                {
                    if (j >= bin_end)
                    {
                        continue;
                    }
                    if (j == pid)
                    {
                        continue;
                    }
                    double rxd = rxbg[j] - rxl;
                    double ryd = rybg[j] - ryl;
                    double rzd = rzbg[j] - rzl;
                    double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
                    double rinv = rsqrt(d2); //distance
                    double d = rinv*d2;
                    //kernel returns pair. first element is energy, 2nd is dvdr/r
                    thrust::pair<double, double> energy_dvdroverR = NKernel2(p, d, rinv, rinv * rinv);

                    //multiply energy by .5 because each pair is process twice in this alg.
                    ee += thrust::get<0>(energy_dvdroverR);
                    ffx += thrust::get<1>(energy_dvdroverR) * rxd;
                    ffy += thrust::get<1>(energy_dvdroverR) * ryd;
                    ffz += thrust::get<1>(energy_dvdroverR) * rzd;
                    //printf("fx %f\n", ffx); 
                    //printf("i%i j%i e %f   \n",pid,j,.5*thrust::get<0>(energy_dvdroverR) );
                    CHECK(
                          if (pid == check)
                    {
                          //a printf("e %f   \n",thrust::get<0>(energy_dvdroverR) );
                          printf("j! %i rx %f ry %f\n", j, rxbg[j], rybg[j]);
                    }
                          )

                }
                CHECK(if (pid == check)printf("\n");)
                    //}//
                    //todo check if the *2 is ok
                }
        }
    }
    fx[ii] += ffx;
    fy[ii] += ffy;
    fz[ii] += ffz;
    e[ii] += .5 * ee;
}


//does the force calulation without a neighbor list
//explore all 27 bins and zeros out the forces that are outside of rcut
//currently doesn't use periodic boundary conditions, so only work with open conditions

template <typename K>
__global__ void processPairNaive(int *rback, double *rxbg, double *rybg, double *rzbg,
                                 double *e, double *fx, double *fy, double *fz,
                                 unsigned *bin, int *nbins, int *nbrIdsx, int *nbrIdsy, int *nbrIdsz,
                                 int *binHeads, int* numNbrs, int *binCounts, K * parms, int n)
{

    __shared__ double shm_e[SHARED_EVAL * TPP_EVALp]; //TODO use synch shuffle for e&f to reserve shmem for page tables
    __shared__ double shm_fx[SHARED_EVAL * TPP_EVALp]; //TODO use synch shuffle for e&f to reserve shmem for page tables
    __shared__ double shm_fy[SHARED_EVAL * TPP_EVALp]; //TODO use synch shuffle for e&f to reserve shmem for page tables
    __shared__ double shm_fz[SHARED_EVAL * TPP_EVALp]; //TODO use synch shuffle for e&f to reserve shmem for page tables

    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    int tid = threadIdx.y;
    int tid_x = threadIdx.x;
    //shouldn't have more threads than particles
    int nion = n;
    if (pid >= nion)
    {
        return;
    }

    int ii = rback[pid];
    if (ii >= n) return; //exit it not local
    //int binsTot = nbins[0]*nbins[1]*nbins[2];
    int local_bin = bin[pid];

    //get local particle assigned to this thread
    double rxl = rxbg[pid];
    double ryl = rybg[pid];
    double rzl = rzbg[pid];

    int loc_bin = local_bin;
    //if bin is on the boundary, some neighboring bins wont exist. skip those bins
    //this shouldn't result in branch diversion if each warp hits same condition
    int xybins = nbins[0] * nbins[1];
    int loc_bin_z = local_bin / (xybins);
    int loc_bin_xy = (loc_bin - loc_bin_z * (xybins));
    int loc_bin_y = loc_bin_xy / (nbins[0]);
    int loc_bin_x = loc_bin_xy - loc_bin_y * (nbins[0]);


    double eee = 0;
    double fffx = 0;
    double fffy = 0;
    double fffz = 0;

    //if (pid==0) printf("nx %i ny %i nz %i\n", numNbrs[0], numNbrs[1], numNbrs[2]);
    for (int i = 0; i < numNbrs[0]; i++)
    {
        for (int j = 0; j < numNbrs[1]; j++)
        {
            for (int k = 0; k < numNbrs[2]; k++)
            {

                //if bin is on the boundary, some neighboring bins will be 'out of bounds' (ghost bins)
                //this shouldn't result in branch diversion if each warp hits same condition
                //however, if a bin is out of bounds but also lies in a periodic direction
                //we should process that bin as a ghost bin and shift the particles accordingly
                int nbr_bin_z = loc_bin_z + nbrIdsz[k];
                int nbr_bin_y = loc_bin_y + nbrIdsy[j];
                int nbr_bin_x = loc_bin_x + nbrIdsx[i];


                //if (pid==0) printf("nox %i noy %i noz %i\n", i, j, k);
                //whether neighbor bin is out of bounds, not considering periodicity
                int out_x = nbr_bin_x < 0 || nbr_bin_x >= nbins[0];
                int out_y = nbr_bin_y < 0 || nbr_bin_y >= nbins[1];
                int out_z = nbr_bin_z < 0 || nbr_bin_z >= nbins[2];

                //pbc is commented out for now...
                /*
                //whether neighbor bin is out of bounds, but now considering periodicity
                //ie if bin is out of bounds but in periodic direction (ghost bin), it is valid
                //if not periodic, then the neighbor bin is invalid and not considered later
                //if bin is in bounds, it is trivially valid
                int valid_x =(!out_x) || periodicity[0];
                int valid_y =(!out_y) || periodicity[1];
                int valid_z =(!out_z) || periodicity[2];
      

                //whether bin is a ghost bin
                //a bin is is a ghost bin if it is 'out of bounds' and lies in a periodic direction
                int ghost_x = out_x && periodicity[0];
                int ghost_y = out_y && periodicity[1];
                int ghost_z = out_z && periodicity[2];
     
                //printf("pid %i ghost_x %i periodic x %i\n", pid, ghost_x, periodicity[0]);
                //skip bins that are not in bounds and not in a periodic direction
      
 
                int valid_x =(!out_x);
                int valid_y =(!out_y);
                int valid_z =(!out_z);
                 */

                if (out_x || out_y || out_z)
                {
                    continue;
                }

                /*
                //if bin is a ghost bin, find which bin it is 'replicating'
                if (ghost_x || ghost_y || ghost_z)
                {
                   nbr_bin_x = modulo(nbr_bin_x,nbins[0]);
                   nbr_bin_y = modulo(nbr_bin_y,nbins[1]);
                   nbr_bin_z = modulo(nbr_bin_z,nbins[2]);
                }
                 */
                /*if(pid==0)
                {
                   printf("cur bin %i nbrx %i %i \n",loc_bin,  nbr_bin_x,  nbr_bin_y);
                }*/

                int nbr_bin = nbr_bin_z * nbins[0] * nbins[1] + nbr_bin_y * nbins[0] + nbr_bin_x;
                //now we have to loop over remote particles in neighboring bin.  
                int bin_end = binHeads[nbr_bin] + binCounts[nbr_bin];
                //int chunksize = ceil((float)binCounts[nbr_bin]/TPP_EVAL);
                //int chunk_start = tid_x*chunksize+binHeads[nbr_bin];
                //int chunk_end = min(binHeads[nbr_bin]+binCounts[nbr_bin], (tid_x+1)*chunksize+binHeads[nbr_bin]);
                //binHeads points to start of bin. binCounts is length of bin
                //for (int jj = binHeads[nbr_bin];jj<binHeads[nbr_bin]+binCounts[nbr_bin]; jj++)
                for (int jj = binHeads[nbr_bin] + tid_x; jj < bin_end; jj += TPP_EVAL)
                    //for (int jj =chunk_start;jj<chunk_end; jj+=1)
                {
                    if (jj >= bin_end)
                    {
                        continue;
                    }
                    if (jj == pid)
                    {
                        continue;
                    }
                    double rxd = rxbg[jj] - rxl;
                    double ryd = rybg[jj] - ryl;
                    double rzd = rzbg[jj] - rzl;
                    double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
                    double rinv = rsqrt(d2); //distance
                    double d = rinv*d2;
                    //kernel returns pair. first element is energy, 2nd is dvdr/r
                    thrust::pair<double, double> energy_dvdroverR = NKernel2(*parms, d, rinv, rinv * rinv);
                    //if (pid==0) printf("%i %i r1 %f %f %f r2 %f %f %f\n",pid,jj, rxl, ryl, rzl, rxbg[jj], rybg[jj], rzbg[jj]);
                    //multiply energy by .5 because each pair is process twice in this alg.
                    eee += thrust::get<0>(energy_dvdroverR);
                    fffx += thrust::get<1>(energy_dvdroverR) * rxd;
                    fffy += thrust::get<1>(energy_dvdroverR) * ryd;
                    fffz += thrust::get<1>(energy_dvdroverR) * rzd;

                }
            }
        }
    }
    shm_e[tid * TPP_EVALp + tid_x] = eee;
    shm_fx[tid * TPP_EVALp + tid_x] = fffx;
    shm_fy[tid * TPP_EVALp + tid_x] = fffy;
    shm_fz[tid * TPP_EVALp + tid_x] = fffz;
    __syncthreads();
    double ffx = 0;
    double ffy = 0;
    double ffz = 0;
    double ee = 0;
    /*
       if (tid<16){
          shm_e[tid*TPP_EVALp+tid_x]+=shm_e[tid*TPP_EVALp+tid_x+16];
       }
       __syncthreads();
       if (tid<8){
          shm_e[tid*TPP_EVALp+tid_x]+=shm_e[tid*TPP_EVALp+tid_x+8];
       }
       __syncthreads();
     */
    if (tid_x == 0)
    {
        for (int i = 0; i < TPP_EVAL; i++)
        {
            ee += shm_e[tid * TPP_EVALp + i];
            ffx += shm_fx[tid * TPP_EVALp + i];
            ffy += shm_fy[tid * TPP_EVALp + i];
            ffz += shm_fz[tid * TPP_EVALp + i];
        }
        fx[ii] += ffx;
        fy[ii] += ffy;
        fz[ii] += ffz;
        e[ii] += .5 * ee;
    }
}

//this is a more optimized version of processPairNaive that uses shared memory

template <typename K>
__global__ void processPairShared(int *rback, double *rxbg, double *rybg, double *rzbg,
                                  double *e, double *fx, double *fy, double *fz,
                                  unsigned *bin, int *nbins, int *nbrIdsx, int *nbrIdsy, int *nbrIdsz,
                                  int *binHeads, int* numNbrs, int *binCounts, K * parms, int n)
{
    //int check = 4; 
    //1 thread spawned per particle atm, this could change
    //int pid = blockIdx.x*blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    __shared__ double shm_x[SHARED];
    __shared__ double shm_y[SHARED];
    __shared__ double shm_z[SHARED];

    K p = *parms;
    //int binsTot = nbins[0]*nbins[1]*nbins[2];
    //int local_bin = bin[pid]; 
    int local_bin = blockIdx.x;
    int pid = binHeads[local_bin] + tid;

    CHECK(
          if (pid == check)
    {
          printf("local bin %i\n", local_bin);
    }
          if (pid == check)
    {
          printf("bin_tot %i\n", binsTot);
    }
          )

    //get local particle assigned to this thread
    int bin_end = binHeads[local_bin] + binCounts[local_bin];

    double rxl, ryl, rzl;
    if (pid < bin_end) //verify that pid actually exists
    {
        rxl = rxbg[pid];
        ryl = rybg[pid];
        rzl = rzbg[pid];
    }

    int loc_bin = local_bin;
    int ii;
    if (pid < bin_end)
        ii = rback[pid]; //use back pointer to get original particle index


    //if bin is on the boundary, some neighboring bins wont exist. skip those bins
    //this shouldn't result in branch diversion if each warp hits same condition
    int xybins = nbins[0] * nbins[1];
    int loc_bin_z = local_bin / (xybins);

    int loc_bin_xy = (loc_bin - loc_bin_z * (xybins));

    int loc_bin_y = loc_bin_xy / (nbins[0]);
    int loc_bin_x = loc_bin_xy - loc_bin_y * (nbins[0]);

    CHECK(if (pid == check)
    {
          printf("z y x %i %i %i\n", numNbrs[0], numNbrs[1], numNbrs[2]);
    })

    //loop over 27 possible neighboring bins
    for (int i = 0; i < numNbrs[0]; i++)
    {
        for (int j = 0; j < numNbrs[1]; j++)
        {
            for (int k = 0; k < numNbrs[2]; k++)
            {

                //if bin is on the boundary, some neighboring bins wont exist. skip those bins
                //this shouldn't result in branch diversion if each warp hits same condition
                int nbr_bin_z = loc_bin_z + nbrIdsz[k];
                int nbr_bin_y = loc_bin_y + nbrIdsy[j];
                int nbr_bin_x = loc_bin_x + nbrIdsx[i];
                int nbr_bin = nbr_bin_z * nbins[0] * nbins[1] + nbr_bin_y * nbins[0] + nbr_bin_x;
                //if (pid==0){printf("nbr_bin %i, %i %i %i\n",nbr_bin, nbr_bin_x, nbr_bin_y, nbr_bin_z); }
                if (nbr_bin_x < 0 || nbr_bin_x >= nbins[0] || nbr_bin_y < 0 || nbr_bin_y >= nbins[1] || nbr_bin_z < 0 || nbr_bin_z >= nbins[2])
                {
                    continue;
                }
                CHECK(if (pid == check)
                {
                      printf("g_bin %i, %i %i %i\n", nbr_bin, nbr_bin_x, nbr_bin_y, nbr_bin_z); })
                //now we have to loop over remote particles in neighboring bin.  
                double ee = 0;
                double ffx = 0;
                double ffy = 0;
                double ffz = 0;
                //we process them in "chunks" of size=blocksize     
                int chunks = ceil((float) binCounts[nbr_bin] / blockDim.x);
                //if (pid==0){printf("binC ch %i %i  \n\n", binCounts[nbr_bin], chunks);}
                for (int c = 0; c < chunks; c++)
                {
                    //load stuff into shmem
                    if (c * blockDim.x + tid < binCounts[nbr_bin]) //even if thread isn't assigned to thread (pid nonexistent),it can import shmem
                    {
                        int z = binHeads[nbr_bin] + c * blockDim.x + tid;
                        shm_x[tid] = rxbg[z];
                        shm_y[tid] = rybg[z];
                        shm_z[tid] = rzbg[z];
                    }
                    __syncthreads();

                    //binHeads points to start of bin. binCounts is length of bin
                    int nbr_block_start = binHeads[nbr_bin];
                    //int loc_block_start = binHeads[blockIdx.x];
                    int nbr_bin_end = binHeads[nbr_bin] + binCounts[nbr_bin];
                    for (int j = 0; j < blockDim.x; j++)
                    {
                        int jj = j + nbr_block_start;

                        double rxd = shm_x[j] - rxl; // rxbg[jj]
                        double ryd = shm_y[j] - ryl;
                        double rzd = shm_z[j] - rzl;
                        double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
                        __syncthreads();

                        if (pid >= bin_end)
                        {
                            continue;
                        } //skip if we don't need thread
                        if (jj >= nbr_bin_end)
                        {
                            continue;
                        } //break if we're done with bin
                        if (jj == pid)
                        {
                            continue;
                        } //skip self interactions
                        double d = sqrt(d2); //distance

                        //kernel returns pair. first element is energy, 2nd is dvdr/r
                        thrust::pair<double, double> energy_dvdroverR = Kernel(p, d);

                        //multiply energy by .5 because each pair is process twice in this alg.
                        ee += (.5 * thrust::get<0>(energy_dvdroverR));
                        ffx += thrust::get<1>(energy_dvdroverR) * rxd; //shm_x[j];
                        //printf("ff %f\n", thrust::get<1>(energy_dvdroverR));
                        //printf("rdx %f\n",rxd); 
                        ffy += thrust::get<1>(energy_dvdroverR) * ryd; //shm_y[j];
                        ffz += thrust::get<1>(energy_dvdroverR) * rzd; //shm_z[j];

                        //printf("pid %i j %i rx %f ry %f rz %f sx %f sy %f sz %f e %f ee %f ii %i\n",pid,jj, rxbg[jj], rybg[jj], rzbg[jj],shm_x[j], shm_y[j], shm_z[j],.5*thrust::get<0>(energy_dvdroverR),ee,ii );
                    }
                }
                if (pid < bin_end)
                {
                    fx[ii] += ffx;
                    fy[ii] += ffy;
                    fz[ii] += ffz;
                    e[ii] = e[ii] + ee;
                }
                //printf("eee ii %f %i\n", ee, ii);
            }
        }
    }

}


//I'm sure there's a better way to do this
//determines # particles in every bin

__global__ void getBinCounts(int *binCounts, int*listId, unsigned *bin, int n)
{

    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= n)
    {
        return;
    }
    int binID = bin[pid];

    //store start index of particle relative to start of bin
    listId[pid] = atomicAdd(binCounts + binID, 1);
}



//determine the new index of each particle based on
//the particle's bin id
//based off of Tomas Oppelstrup's bin sort cuda code
//except all global synchs are done with cuda 9 cooperative groups
//ie 1 kernel to rule them all

//STAGE 1: block the bins, find partial sum of particles in each block

__global__ static void generateBinPermutation1(int * binCounts, int *partial_sums, int *nbins)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
    int total_bins = nbins[0] * nbins[1] * nbins[2];
    CHECK(if (tid == 0)
    {
          printf("b %i \n", total_bins);
    })

    __shared__ int shm[SHARED];
    //fix this later with padding or something so we don't have to init share mem
    //this way, if we access rogue shmem we are safe
    shm[tid] = 0;
    __syncthreads();
    if (pid < total_bins)
    {
        shm[tid] = binCounts[pid];
    }
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 32; s = s >> 1)
    {
        if (tid < s)
        {
            shm[tid] += shm[tid + s];
        }
        __syncthreads();
        CHECK(if (tid == 0)
        {
              printf("s %i \n", shm[tid]);
        })
    }
    if (tid < 32)
    {
        shm[tid] += shm[tid + 32];
        __syncthreads();
        shm[tid] += shm[tid + 16];
        __syncthreads();
        shm[tid] += shm[tid + 8];
        __syncthreads();
        shm[tid] += shm[tid + 4];
        __syncthreads();
        shm[tid] += shm[tid + 2];
        __syncthreads();
        shm[tid] += shm[tid + 1];
        __syncthreads();
    }
    if (tid == 0)
    {
        partial_sums[blockIdx.x] = shm[0]; //store partial sum to block
        CHECK(if (tid == 0)
        {
              printf("psum %i \n", shm[0]);
        })
    }
}

//GLOBAL SYNCH
//STAGE 2: Compute scan of block counts.
//this way we get the end idx of each block of bins
//this kernel is run in only 1 thread blocka
// when partial_sums exceed size of SHARED the kernel doesn't work
/*
__global__ static void generateBinPermutation2old(int *partial_sums) {
    int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int j = 1;
    __shared__ int shm [SHARED];
    int x = 0;
    shm[tid] = 0;
    shm[tid] = partial_sums[pid];
    while (j < blockDim.x) {
        CHECK(if (tid == 0) {
            printf("j %i \n", j);
        })
        if (tid >= j) {
            x = shm[tid - j];
        }
        __syncthreads();
        if (tid >= j) {
            shm[tid] += x;
        }
        __syncthreads();
        j *= 2;
    }
    partial_sums[pid] = shm[tid];
}
 */

// Call this routine with blockSize <= SHARED, and 1 block (gridSize = 1).

__global__ static void generateBinPermutation2(int *partial_sums, int veclen)
{
    int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;

    // Code will work correctly with more than one block,
    // but won't be any faster.
    if (pid != tid)
    {
        if (tid == 0)
        {
            static int first = 1;
            if (first)
            {
                printf("Should only call this with one block!\n");
                first = 0;
            }
        }
        return;
    }

    __shared__ int shm [SHARED];

    for (int offset = 0; offset < veclen; offset += blockDim.x)
    {
        int j = 1;
        shm[tid] = partial_sums[offset + tid];
        __syncthreads();
        while (j < blockDim.x)
        {
            int x;
            CHECK(if (tid == 0)
            {
                  printf("j %i \n", j);
            })
            if (tid >= j)
            {
                x = shm[tid - j];
            }
            __syncthreads();
            if (tid >= j)
            {
                shm[tid] += x;
            }
            __syncthreads();
            j *= 2;
        }

        const int carry = (offset > 0) ? partial_sums[offset - 1] : 0;
        partial_sums[offset + tid] = shm[tid] + carry;
        __syncthreads();
    }
}

//GLOBAL SYNCH
//STAGE 3: do inclusive  prefix scan within each block, add #items in prev blocks inclusive, and subtract
//#items in current block and #items in current bin to do a bin wide exclusive prefix scan
//after this stage, count array holds an index to where this bin begins

__global__ static void generateBinPermutation3(int *binCounts, int * partial_sums, int *nbins)
{
    int totalbins = nbins[0] * nbins[1] * nbins[2];
    int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ volatile struct
    {
        int psum, x[2][SHARED], y[2][SHARED];
    } shm;

    if (tid == 0) shm.psum = partial_sums[blockIdx.x];

    shm.x[0][tid] = 0;
    shm.y[0][tid] = 0;
    __syncthreads();
    int x = 0;
    if (pid < totalbins)x = binCounts[pid];
    int c = x;
    shm.x[1][tid] = x;
    __syncthreads();

    //prefix scan
    x += shm.x[1][tid - 1];
    shm.y[1][tid] = x;
    __syncthreads();
    x += shm.y[1][tid - 2];
    shm.x[1][tid] = x;
    __syncthreads();
    x += shm.x[1][tid - 4];
    shm.y[1][tid] = x;
    __syncthreads();
    x += shm.y[1][tid - 8];
    shm.x[1][tid] = x;
    __syncthreads();
    x += shm.x[1][tid - 16];
    shm.y[1][tid] = x;
    __syncthreads();
    x += shm.y[1][tid - 32];
    shm.x[1][tid] = x;
    __syncthreads();
    x += shm.x[1][tid - 64];
    shm.y[1][tid] = x;
    __syncthreads();

    //subtractions to make prefix scan exclusive
    CHECK(if (tid == 0)
    {
          printf("pre %i\n", shm.y[1][blockDim.x - 1]);
    })
    x += (shm.psum - shm.y[1][blockDim.x - 1]);
    binCounts[pid] = x - c;
}

//STAGE 4: generate the new indices

__global__ static void generateBinPermutation4(int *binCounts, int *listIds, unsigned* bins, int n)
{
    //int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= n)
    {
        return;
    }
    int bin = bins[pid];
    listIds[pid] += binCounts[bin];
}

void binSort(SYSTEM *sys, int totalbins) 
{
    ACCELERATOR *accelerator = accelerator_getAccelerator(NULL);
    GPUCUDAPARMS *accParms= (GPUCUDAPARMS*)accelerator->parms; 
    int nIon = sys->nion;
    GPUNLIST *gnlist=sys->collection->gnlist;
    int *nbins = gnlist->nbinsg;
    int *binCounts = gnlist->binCountsg;
    int *partial_sums = gnlist->partial_sumsg;
    int *listIds = gnlist->listIdsg;
    unsigned *bins = gnlist->binsg;
    //for kernels iterating over particles
    int blocksize_p = SHARED;
    int gridsize_p = ceil((float) nIon / blocksize_p);

    //for kernels iterating over bins
    int blocksize_b = SHARED;
    int gridsize_b = ceil((float) totalbins / blocksize_b);
    //printf("tot bins %i\n", totalbins);
    //printf("blocksize b %i g size b %i\n", blocksize_b, gridsize_b);
    //for kernels iterating over blocks
    int blocksize_1 = SHARED;
    int gridsize_1 = 1;

    //reduce bins within each block
    //bin counts psums nbins
    generateBinPermutation1 << <gridsize_b, blocksize_b>>>(binCounts, partial_sums, nbins);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    //prefix scan block countsa
    CHECK
        (
         int * psums = (int *) malloc(sizeof (int)*gridsize_b);
         gpu_memcpy_device2host(psums, partial_sums, gridsize_b);
         for (int i = 0; i < gridsize_b; i++)
    {
         printf("psum %i %i \n", i, psums[i]);
    }
         )


    CHECK(
          int nion = nIon;
          int *r_back = (int*) malloc(nion * sizeof (int));
          gpu_memcpy_device2host(r_back, gnlist->listIdsg, nion);
          for (int i = 0; i < nIon; i++)
    {
          printf("%i listid %i\n", i, r_back[i]);

    }
          )

    //printf("perm2 \n");
    
    generateBinPermutation2<<<gridsize_1, blocksize_1>>>(partial_sums, accParms->partial_sumSize);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    CHECK(
          gpu_memcpy_device2host(psums, partial_sums, gridsize_b);
          for (int i = 0; i < gridsize_b; i++)
    {
          printf("psum pre %i %i \n", i, psums[i]);
    }
          )

    //printf("perm3 \n");
    //prefix scan within each block
    generateBinPermutation3 << <gridsize_b, blocksize_b>>>(binCounts, partial_sums, nbins);
    CUDA_SAFE_CALL(cudaPeekAtLastError());

    CHECK(
          int * binScans = (int*) malloc(sizeof (int)*totalbins);
          gpu_memcpy_device2host(binScans, binCounts, totalbins);
          for (int i = 0; i < totalbins; i++)
    {
          printf("scan block %i  %i s %i\n", i, i / blocksize_b, binScans[i]);
    }
          )

    generateBinPermutation4 << <gridsize_p, blocksize_p>>>(binCounts, listIds, bins, nIon);
    CUDA_SAFE_CALL(cudaPeekAtLastError());

    CHECK(
          gpu_memcpy_device2host(r_back, gnlist->listIdsg, nIon);
          for (int i = 0; i < nIon; i++)
    {
          printf("%i newid %i \n", i, r_back[i]);

    }
          )


}

__global__ void calcBinSize(double * mmgpu, double *mins, double *lens, int *nbins, double cutoff) 
{
    int tid = threadIdx.x;
    mins[tid] = mmgpu[2 * tid];
    lens[tid] = mmgpu[2 * tid + 1] - mmgpu[2 * tid];
    nbins[tid] = max(1, (int) floor(lens[tid] / (cutoff)));
}

void binParticlesGPU(SYSTEM * sys, double rcut) 
{

    CUDA_SAFE_CALL(cudaPeekAtLastError());
    GPUNLIST *gnlist=sys->collection->gnlist;
    //int nLocal = sys->nlocal; 
    int nIon = sys->nion;
    //size_t bytes = nIon*sizeof(double);
    //STATE *state = sys->collection->state;      
    STATE *gsh = sys->collection->gpustate_h;
    //double * rxbg,*rybg, *rzbg;
    int * r_backg;
    unsigned *binsg;
    double *mmgpu = gnlist->mmgpu;
    //this can be put in constant memory
    //rxbg = sys->collection->rxbg; rybg = sys->collection->rybg; rzbg = sys->collection->rzbg;
    r_backg = gnlist->r_backg;
    binsg = gnlist->binsg;

    double mmgpuH[6];
    //double mmgpuH1[6];

    //fill back pointer array
    int blocksize = SHARED;
    int gridsize = ceil((float) nIon / blocksize);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    fillBackPointers << <gridsize, blocksize>>>(r_backg, nIon);
    CUDA_SAFE_CALL(cudaPeekAtLastError());

    //get minimum and  maximum x,y,z coorindates of particles in state
    //if using periodic boundary conditions, use the predefined computational box
    //otherwise, find the min/max coordinates
    if (0 && sys->box->pbc)
    {
        THREE_VECTOR corner = box_get_corner(sys->box);
        THREE_VECTOR corner2 = box_get_urcorner(sys->box);
        mmgpuH[0] = corner.x;
        mmgpuH[1] = corner2.x;
        mmgpuH[2] = corner.y;
        mmgpuH[3] = corner2.y;
        mmgpuH[4] = corner.z;
        mmgpuH[5] = corner2.z;
        gpu_memcpy_host2device(mmgpu, mmgpuH, 6);
    }
    else
    {
        getMinMaxGPU(sys->collection, sys->collection->gpustate_h, sys->collection->state, mmgpu);
        gpu_memcpy_device2host(mmgpuH, mmgpu, 6);
        for (int i = 0; i < 6; i++)
        {
            CHECK_1(printf("Minmax %i %f\n", i, mmgpuH[i]);)
                //mmgpuH1[i] = mmgpuH[i];
        }
    }

    //double skin = sys->neighbor->deltaR;


    calcBinSize << <1, 3 >> >(mmgpu, gnlist->minsg, gnlist->lensg, gnlist->nbinsg, rcut);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    //determine size of sim box
    double lens[3];
    int nbins[3];
    //int mins[3];
    //double binL[3];

    lens[0] = mmgpuH[1] - mmgpuH[0];
    lens[1] = mmgpuH[3] - mmgpuH[2];
    lens[2] = mmgpuH[5] - mmgpuH[4];
    //mins[0] = mmgpuH[0];
    //mins[1] = mmgpuH[2];
    //mins[2] = mmgpuH[4];

    double *minsg, *lensg;
    minsg = gnlist->minsg;
    lensg = gnlist->lensg;
  
    //determine number of bins necessary based on cutoff
    CHECK_1(
            printf("box len %f %f %f \n", lens[0], lens[1], lens[2]);
            )

        //skin = sys->neighbor->deltaR; //.01; 
        nbins[0] = std::max(1, (int) floor(lens[0] / (rcut)));
    nbins[1] = std::max(1, (int) floor(lens[1] / (rcut)));
    nbins[2] = std::max(1, (int) floor(lens[2] / (rcut)));
    //binL[0] = lens[0] / nbins[0];
    //binL[1] = lens[1] / nbins[1];
    //binL[2] = lens[2] / nbins[2];
    int nBinsTot = nbins[0] * nbins[1] * nbins[2];
    gnlist->nBinsTot = nBinsTot;
    //int tBins = nbins[0]*nbins[1]*nbins[2]; 

    blocksize = SHARED;
    gridsize = ceil((float) nIon / blocksize);

    //Chceck the change of total number of bins and re-alloc GPU memory accordingly.
    allocGPUnbins(gnlist, nBinsTot, blocksize);

    CHECK_1(
            printf("nBins %i %i %i \n", nbins[0], nbins[1], nbins[2]);
            //printf("binLens %f %f %f\n", binL[0], binL[1], binL[2]);
            printf("gridsize %i blocksize %i nion %i\n", blocksize, gridsize, nIon);
            )


        //determine binid for each particles
        assignBinIDs << <gridsize, blocksize>>>(binsg, gsh->rx, gsh->ry, gsh->rz, lensg, gnlist->nbinsg, minsg, nIon);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    CHECK_3(
            unsigned *bin_x = (unsigned*) malloc(nIon * sizeof (unsigned));
            gpu_memcpy_device2host(bin_x, binsg, nIon);
            for (int i = 0; i < nIon; i++)
    {
            printf("i %i %f %f %f bin %i \n", i, state->rx[i], state->ry[i], state->rz[i], bin_x[i]);
    }
            free(bin_x);
            )

        //count particles/bin
        CUDA_SAFE_CALL(cudaMemset(gnlist->binCountsg, 0, nBinsTot * sizeof (int)));
    CUDA_SAFE_CALL(cudaMemset(gnlist->listIdsg, 0, nIon * sizeof (int)));
    getBinCounts << <gridsize, blocksize>>>(gnlist->binCountsg, gnlist->listIdsg, binsg, nIon);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    gpu_memcpy_device2device(gnlist->binCountsg2, gnlist->binCountsg, nbins[0] * nbins[1] * nbins[2]);

    CHECK_2(
            int *counts = (int*) malloc(nBinsTot * sizeof (int));
            gpu_memcpy_device2host(counts, gnlist->binCountsg, nBinsTot);
            for (int i = 0; i < nbins[0] * nbins[1] * nbins[2]; i++)
    {
            printf("bin %i c %i\n", i, counts[i]);
    }
            free(counts);
            )

        //do a sort on binds to determine new indices of each particle   
        binSort(sys, nBinsTot);
    CUDA_SAFE_CALL(cudaPeekAtLastError());

    //binCountsg is actuall bin heads now
    gnlist->binHeadsg = gnlist->binCountsg;

    //bin particles by assigning particles to new indices
    cudaMemset(gnlist->binsbg, 0, nIon * sizeof (unsigned));
    permuteParticles2(sys->collection, nIon); //wraps the kernel that populates the bin-sorted arrays of particles
    CUDA_SAFE_CALL(cudaPeekAtLastError());

    CHECK(
          unsigned *bin_x = (unsigned*) malloc(sizeof (unsigned)*nIon);
          gpu_memcpy_device2host(bin_x, gnlist->binsbg, nIon);
          for (int i = 0; i < nIon; i++)
    {
          printf("bin %i\n", bin_x[i]);
    }
          int * r_back1 = (int *) malloc(nIon * sizeof (int));
          int * r_backbg = (int *) malloc(nIon * sizeof (int));
          gpu_memcpy_device2host(r_back1, gnlist->listIdsg, nIon);
          gpu_memcpy_device2host(r_backbg, gnlist->r_backg, nIon);
          for (int i = 0; i < nIon; i++)
    {
          printf("%i r %i rb %i\n", i, r_back1[i], r_backbg[i]);
    }
          free(bin_x);
          )

        CUDA_SAFE_CALL(cudaPeekAtLastError());
    CHECK_2(
            counts = (int *) malloc(sizeof (int)*nBinsTot);
            gpu_memcpy_device2host(counts, (gnlist->binCountsg2), nBinsTot);
            int *heads = (int *) malloc(sizeof (int)*nBinsTot);
            gpu_memcpy_device2host(heads, (gnlist->binHeadsg), nBinsTot);
            for (int i = 0; i < nBinsTot; i++)
    {
            printf("binC %i %i\n", i, counts[i]);
    }
            for (int i = 0; i < nBinsTot; i++)
    {
            printf("binH %i %i\n", i, heads[i]);
    }
            free(counts);
            )



        //determine offsets for determining 27 potential neighbors of a bin
        //27 if 1d
        // 9 if 2d
        // 3 if 1d
        int* x_nbr=NULL;
    int* y_nbr=NULL;
    int* z_nbr=NULL;
    int num_x_nbr = 3;
    int num_y_nbr = 3;
    int num_z_nbr = 3;
    if (nbins[0] > 1)
    {
        //int x_nbr1 [3] = {-1, 0, 1};
        //x_nbr = x_nbr1;
        x_nbr=(int *)ddcMalloc(num_x_nbr*sizeof(int));
        x_nbr[0]=-1;
        x_nbr[1]=0;
        x_nbr[2]=1;
    }
    else
    {
        //int x_nbr1 [1] = {0};
        //x_nbr = x_nbr1;
        num_x_nbr = 1;
        x_nbr=(int *)ddcMalloc(num_x_nbr*sizeof(int));
        x_nbr[0]=0;
    }
    if (nbins[1] > 1)
    {
        //int y_nbr1 [3] = {-1, 0, 1}; //{-nbins[0],0,nbins[0]};
        //y_nbr = y_nbr1;
        y_nbr=(int *)ddcMalloc(num_y_nbr*sizeof(int));
        y_nbr[0]=-1;
        y_nbr[1]=0;
        y_nbr[2]=1;

    }
    else
    {
        //int y_nbr1 [1] = {0};
        //y_nbr = y_nbr1;
        num_y_nbr = 1;
        y_nbr=(int *)ddcMalloc(num_y_nbr*sizeof(int));
        y_nbr[0]=0;
    }
    if (nbins[2] > 1)
    {
        //int z_nbr1 [3] = {-1, 0, 1}; //{-nbins[0]*nbins[1],0,nbins[0]*nbins[1]};
        //z_nbr = z_nbr1;
        z_nbr=(int *)ddcMalloc(num_z_nbr*sizeof(int));
        z_nbr[0]=-1;
        z_nbr[1]=0;
        z_nbr[2]=1;
    }
    else
    {
        //int z_nbr1 [1] = {0};
        //z_nbr = z_nbr1;
        num_z_nbr = 1;
        z_nbr=(int *)ddcMalloc(num_z_nbr*sizeof(int));
        z_nbr[0]=0;
    }

    CHECK_1(
            for (int i = 0; i < num_x_nbr; i++)
    {
            printf("x nbrs %i\n", x_nbr[i]);
    }
            for (int i = 0; i < num_y_nbr; i++)
    {
            printf("y nbrs %i\n", y_nbr[i]);
    }
            for (int i = 0; i < num_z_nbr; i++)
    {
            printf("z nbrs %i\n", z_nbr[i]);
    }
            )

    int num_nbrsxyz[3];
    num_nbrsxyz[0] = num_x_nbr;
    num_nbrsxyz[1] = num_y_nbr;
    num_nbrsxyz[2] = num_z_nbr;

    int num_nbrs = num_z_nbr * num_y_nbr*num_x_nbr;
    gnlist->numNbrs = num_nbrs;

    if (0)
    {
        int nbrs[num_nbrs]; // = (int*) malloc(num_nbrs*sizeof(int));
        int nbrsxy = num_x_nbr*num_y_nbr;
        for (int i = 0; i < num_z_nbr; i++)
        {
            for (int j = 0; j < num_y_nbr; j++)
            {
                for (int k = 0; k < num_x_nbr; k++)
                {
                    nbrs[num_x_nbr * num_y_nbr * i + num_x_nbr * j + k] = z_nbr[i] + y_nbr[j] + x_nbr[k];
                    printf("nbrs %i %i %i %i\n", z_nbr[i], y_nbr[j], x_nbr[k], nbrs[nbrsxy * i + num_x_nbr * j + k]);
                }
            }
        }
    }

    gpu_memcpy_host2device(gnlist->nbrIdsx, x_nbr, num_x_nbr);
    gpu_memcpy_host2device(gnlist->nbrIdsy, y_nbr, num_y_nbr);
    gpu_memcpy_host2device(gnlist->nbrIdsz, z_nbr, num_z_nbr);
    gpu_memcpy_host2device(gnlist->numNbrsxyz, num_nbrsxyz, 3);


    // Free the memory of pointers
    ddcFree(x_nbr);
    ddcFree(y_nbr);
    ddcFree(z_nbr);
}

//
//gpu version of non neighbor list lennard jones pair processor that uses bins
//follows same interface as other ddcmd potential functions

void pairEvaluateNaiveGPU(SYSTEM *sys, PAIR_PARMS *parms, ETYPE *e) 
{
    GPUNLIST *gnlist=sys->collection->gnlist;
    STATE *gsh = sys->collection->gpustate_h;
    LJGPU *lj = (LJGPU*) malloc(sizeof (LJGPU));

    LJ_PARMS** ljp = (LJ_PARMS**) parms->parms;
    lj->eps = ljp[0]->eps;
    lj->sig = ljp[0]->sigma;
    lj->shift = ljp[0]->shift;
    lj->rcut = parms->rcut[0].value;

    LJGPU *ljgpu;
    gpu_allocator(ljgpu, 1);
    gpu_memcpy_host2device(ljgpu, lj, 1);
    int nion = sys->nion;
    double *e_all = gnlist->e_all;
    cudaMemset(e_all, 0, sizeof (double)*sys->nion);
    cudaMemset(gsh->fx, 0, sizeof (double)*sys->nion);
    cudaMemset(gsh->fy, 0, sizeof (double)*sys->nion);
    cudaMemset(gsh->fz, 0, sizeof (double)*sys->nion);

    //launch neighbor list force/energy kernel
    int blockSize = SHARED_EVAL;

    int gridSize = ceil((float) nion / blockSize);
    dim3 blockSize3(TPP_EVAL, blockSize);

    CHECK(printf("launching process pair naive\n");)
    processPairNaive << <gridSize, blockSize3>>>(gnlist->r_backg, gnlist->rxbg, gnlist->rybg, gnlist->rzbg, e_all,
        gsh->fx, gsh->fy, gsh->fz, gnlist->binsbg,
        gnlist->nbinsg, gnlist->nbrIdsx, gnlist->nbrIdsy, gnlist->nbrIdsz,
        gnlist->binHeadsg, gnlist->numNbrsxyz, gnlist->binCountsg2, ljgpu, sys->nion);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
}

void pairEvaluateSharedGPU(SYSTEM *sys, PAIR_PARMS *parms, ETYPE *e) 
{
    STATE *gsh = sys->collection->gpustate_h;
    GPUNLIST *gnlist=sys->collection->gnlist;

    //send stuff to gpu 
    LJGPU *lj = (LJGPU*) malloc(sizeof (LJGPU));
    LJ_PARMS** ljp = (LJ_PARMS**) parms->parms;
    lj->eps = ljp[0]->eps;
    lj->sig = ljp[0]->sigma;
    lj->shift = ljp[0]->shift;
    lj->rcut = parms->rcut[0].value;
    LJGPU *ljgpu;
    gpu_allocator(ljgpu, 1);
    gpu_memcpy_host2device(ljgpu, lj, 1);

    //zero energy buffer   
    double *e_all = gnlist->e_all;
    CUDA_SAFE_CALL(cudaMemset(gnlist->e_all, 0, sizeof (double)*sys->nlocal));

    //launch kernel
    int blocksize = SHARED;
    int gridsize = gnlist->nBinsTot;
    CHECK(printf("launching process pair naive\n");)
    processPairShared << <gridsize, blocksize>>>(gnlist->r_backg, gnlist->rxbg, gnlist->rybg, gnlist->rzbg, e_all,
        gsh->fx, gsh->fy, gsh->fz, gnlist->binsbg,
        gnlist->nbinsg, gnlist->nbrIdsx, gnlist->nbrIdsy, gnlist->nbrIdsz,
        gnlist->binHeadsg, gnlist->numNbrsxyz, gnlist->binCountsg2, ljgpu, sys->nion);
    CUDA_SAFE_CALL(cudaPeekAtLastError());

}

#ifdef __cplusplus
extern "C"
{
#endif

#define GPU_REDUCE 1
#define CPU_REDUCE 0
#define DEBUG_REDUCE(x) 

void calcTensorWithGPUVirials(SYSTEM *sys, ETYPE *e) 
{
    int nIon = sys->nion;
    int nLocal = sys->nlocal;
    GPUNLIST *gnlist=sys->collection->gnlist;
    double *sionxx = gnlist->sion_h->xx;
    double *sionyy = gnlist->sion_h->yy;
    double *sionzz = gnlist->sion_h->zz;
    //double *sionxy = gnlist->sion_h->xy;
    //double *sionxz = gnlist->sion_h->xz;
    //double *sionyz = gnlist->sion_h->yz;
    double virSums[3];
    double virCors[3];
    virSums[0] = 0;
    virSums[1] = 0;
    virSums[2] = 0;

    if (CPU_REDUCE)
    {
        gpu_memcpy_device2host(sionxx, gnlist->gpu_sion->xx, nIon);
        gpu_memcpy_device2host(sionyy, gnlist->gpu_sion->yy, nIon);
        gpu_memcpy_device2host(sionzz, gnlist->gpu_sion->zz, nIon);

        double ssxx = 0;
        double ssyy = 0;
        double sszz = 0; // double ssxy=0; double ssxz=0; double ssyz=0;
        for (int i = 0; i < sys->nion; i++)
        {
            ssxx -= sionxx[i];
            ssyy -= sionyy[i];
            sszz -= sionzz[i];
        }
        DEBUG_REDUCE(printf("gpu->cpu sions %f %f %f\n", ssxx, ssyy, sszz);)
        virSums[0] = ssxx;
        virSums[1] = ssyy;
        virSums[2] = sszz;
    }
    if (GPU_REDUCE)
    {
        double* virSums_d = NULL;
        gpu_allocator(virSums_d, 3);
        //size_t avail_bytes = gnlist->allocFactor*nLocal*sizeof(double);
        size_t temp_storage_bytes = 0;
        void *d_temp_storage = NULL;
        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->gpu_sion->xx, virSums_d, nLocal);
        cudaMalloc(&d_temp_storage, temp_storage_bytes);

        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->gpu_sion->xx, virSums_d, nLocal);
        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->gpu_sion->yy, virSums_d + 1, nLocal);
        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->gpu_sion->zz, virSums_d + 2, nLocal);
        double* virSumsPtr = virSums;
        gpu_memcpy_device2host(virSumsPtr, virSums_d, 3);
        virSums[0] *= -1;
        virSums[1] *= -1;
        virSums[2] *= -1;
        DEBUG_REDUCE(printf("gpu sions %8.15f %8.15f %8.15f\n", virSums[0], virSums[1], virSums[2]);)

        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->virCorx, virSums_d, nLocal);
        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->virCory, virSums_d + 1, nLocal);
        cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->virCorz, virSums_d + 2, nLocal);
        double* virCorsPtr = virCors;
        gpu_memcpy_device2host(virCorsPtr, virSums_d, 3);
        virSums[0] -= virCors[0];
        virSums[1] -= virCors[1];
        virSums[2] -= virCors[2];
        DEBUG_REDUCE(printf("c gpu sions %8.15f %8.15f %8.15f\n", virSums[0], virSums[1], virSums[2]);)

        if (d_temp_storage) cudaFree(d_temp_storage);
        if (virSums_d) cudaFree(virSums_d);

        CUDA_SAFE_CALL(cudaPeekAtLastError());


    }

    e->sion.xx = virSums[0];
    e->sion.yy = virSums[1];
    e->sion.zz = virSums[2];

    e->virial.xx = virSums[0];
    e->virial.yy = virSums[1];
    e->virial.zz = virSums[2];
}
#undef GPU_REDUCE 
#undef CPU_REDUCE 
#undef DEBUG_REDUCE

__global__ void moveEnergies(double *e, double *results, int nLocal)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nLocal)return;
    e[pid] = results[7 * pid];
}

void sendForceEnergyToHost(SYSTEM *sys, ETYPE *e) 
{
    STATE * state = sys->collection->state;
    //STATE *gsh = sys->collection->gpustate_h;
    GPUNLIST *gnlist=sys->collection->gnlist;
    //int nIon = state->nion;
    int nLocal = state->nion;

    int blockSize = 32;
    int gridSize = ceil((float) nLocal / blockSize);
    gridSize = ceil((float) sys->nlocal / blockSize);
    //  moveEnergies<<<gridSize, blockSize>>>(sys->collection->e_all, gnlist->results, sys->nlocal);
    CUDA_SAFE_CALL(cudaPeekAtLastError(););
    //cudaThreadSynchronize();
    cudaDeviceSynchronize();
    //double esum=0;
    //double esum1=0;
    double esumAll[2];

    double* esumAll_d = NULL;
    gpu_allocator(esumAll_d, 2);

    //reduce energy
    size_t temp_storage_bytes = 0;
    void *d_temp_storage = NULL;
    //size_t avail_bytes = gnlist->allocFactor*nLocal*sizeof(double);
    //int DSIZE=sys->collection->state->nion*sizeof(double);
    //double *h_data = (double *)malloc(DSIZE);
    //cudaMemcpy(h_data, sys->collection->e_all, DSIZE, cudaMemcpyDeviceToHost);
    //for(int i=0; i<sys->collection->state->nion; i++){
    //if(h_data[i]>0.1) printf("i=%d   e_all=%f\n", i, h_data[i]);
    //}
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
    cudaMalloc(&d_temp_storage, temp_storage_bytes);

    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e_all, esumAll_d, nLocal);
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gnlist->e1, esumAll_d + 1, nLocal);

    gpu_memcpy_device2host(esumAll, esumAll_d, 2);
    e->eion += (esumAll[0] + esumAll[1]);

    if (d_temp_storage) cudaFree(d_temp_storage);
    if (esumAll_d) cudaFree(esumAll_d);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    //calcMolecularPressuresGPU(sys);
    //calcTensorWithGPUVirials(sys, e);
}

void sendEnergyToHost(SYSTEM *sys, ETYPE *e) 
{
    GPUNLIST *gnlist = sys->collection->gnlist;
    //STATE *gsh = sys->collection->gpustate_h;
    STATE * state = sys->collection->state;
    int nIon = state->nion;
    double *ecpu = (double *) malloc(sizeof (double)*sys->nion);
    memset(ecpu, 0, sizeof (double)*sys->nion);
    e->eion = 0;
    //printf("nion %i\n", nIon);   
    gpu_memcpy_device2host(ecpu, gnlist->e_all, nIon);

    //TODO this reduction can also be done on GPU
    //trivial, will get back to this
    double ee = 0;
    for (int i = 0; i < sys->nion; i++)
    {
        //printf("i %i e %f\n", i, ecpu[i]);
        ee += ecpu[i];
    }
    e->eion += ee;
    //printf("eion %f\n", e->eion);
}

//bins particles, and evaluates forces

void pairProcessTemplatedGpu(SYSTEM *sys, PAIR_PARMS *parms, ETYPE *e) 
{
    //struct timeval start_send, end_send, start_bin, end_bin, start_pair, end_pair; 

    //GPROF(gettimeofday(&start_send, NULL);)
    //send state from host to device
    sendGPUState(sys, sys->nlocal);
    //GPROF(cudaStreamSynchronize(0);
    //gettimeofday(&end_send, NULL);)

    //GPROF(gettimeofday(&start_bin, NULL);)
    double rcut = parms->rcut[0].value;
    //COLLECTION * collection = sys->collection;
    //int nion = sys->nion;
    //sort particles into bins
    binParticlesGPU(sys, rcut);
    //GPROF(cudaStreamSynchronize(0);
    //gettimeofday(&end_bin, NULL);)

    //GPROF(gettimeofday(&start_pair, NULL);)
    //evaluate all possible pairs using shared memory kernel
    //pairEvaluateSharedGPU(sys, parms, e); //uses shared memory
    //permuteParticles2(collection, nion);
    pairEvaluateNaiveGPU(sys, parms, e); //doesn't use shared memory
    //GPROF(cudaStreamSynchronize(0);
    //gettimeofday(&end_pair, NULL);)
    /*
    GPROF(
    printf("gpu send time %ld\n", ((end_send.tv_sec * 1000000 + end_send.tv_usec)
                   - (start_send.tv_sec * 1000000 + start_send.tv_usec)));

    printf("gpu bin %ld\n", ((end_bin.tv_sec * 1000000 + end_bin.tv_usec)
                   - (start_bin.tv_sec * 1000000 + start_bin.tv_usec)));

    printf("gpu pair time %ld\n", ((end_pair.tv_sec * 1000000 + end_pair.tv_usec)
                   - (start_pair.tv_sec * 1000000 + start_pair.tv_usec)));
         )
     */
    //send forces and energy from gpu to host
    sendForceEnergyToHost(sys, e);

}
#ifdef __cplusplus
}
#endif


