#include <cassert>
#include "nlistGPU.h"
#include "cudaUtils.h"
#include <cuda_fp16.h>
#include <thrust/pair.h>
#include "pairProcessGPU.h"
#include "bioCharmmParms.h"
#include "bioMartini.h"
#include "units.h"
#include "codata.h"
//#include "heap.h"
#include "ddcMalloc.h"
#include "gpuMemUtils.h"
#include "cudaTypes.h"
#include "simulate.h"
#include "gpu_allocator.hpp"
#include "accelerator.h"

//debugging only
//#include "cuda_profiler_api.h"
#define USE_DVIRIAL
#define CHECK(x)
#ifdef USE_DVIRIAL
#define DVIRIAL(x) x
#define NUMBUFFS 7
#else
#define NUMBUFFS 4
#define DVIRIAL(x) 
#endif
#define TPP_NLIST_EVAL 8//threads per particle used in nlist evaluator kernel
#define TPP_NLIST_CONSTRUCT 4
#define SHARED_NLIST 4  //shared mem perblock in nlist kernels
#define PBC(x) x   
#define VIRIAL(x)   
//useful for debuging, not used otherwise:
//static int nnpages;
static int nnppgpu;
//static int ppgpu;
//const static int maxPagesPerParticle = 1;

/*
this function allocates a pool of pages of gpu memory
this way, when building a neibhor list
we can quickly resize each particle's neighbor list buffer
by adding or removing pre allocated to its list when necessary
TODO: params like page size, expected max num of pages per particle, 
and expected max #particles/gpu should be adjustable in object.data file
 */
void allocPages(GPUNLIST * gpunlist, int particlesPerGpu, int density, double cutoff) 
{
    ACCELERATOR *accelerator = accelerator_getAccelerator(NULL);
    GPUCUDAPARMS *accParms = (GPUCUDAPARMS *) accelerator->parms; 
    printf("particles per gpu %i\n", particlesPerGpu);
    //again, we allocate 1.5 times as much memory to be safe
    //particlesPerGpu*=1.2;
    //	printf("particles per gpu up %i\n", particlesPerGpu);
    //ppgpu = particlesPerGpu;
    //GPUNLIST * gnl = gpunlist;
    //double PI=3.141592654;
    //double volume = (4/3)*PI*cutoff*cutoff*cutoff; //vol of cutoff sphere
    //int maxPagesPerParticle = 16; 
    int npages = ceil((float) particlesPerGpu / PAGE_SIZE);

   
    gpunlist->plist_h = (int**) malloc(sizeof (int*)*npages); //alloc gpu page list on cpu
    gpu_allocator(gpunlist->plist, npages);
    gpu_allocator(gpunlist->nbrSizes, particlesPerGpu);
    gpu_allocator(gpunlist->mem, accParms->maxPagesPerParticle * PAGE_SIZE * particlesPerGpu);

    int numResultBuffers = 4;
    DVIRIAL(numResultBuffers = 7);
    gpu_allocator(gpunlist->results, numResultBuffers * particlesPerGpu);

    //used for debugging
    //nnpages = npages;
    nnppgpu = particlesPerGpu;

    //create pointers to each page
    for (int i = 0; i < npages; i++) {
        gpunlist->plist_h[i] = gpunlist->mem + PAGE_SIZE*i;
    }

    //send array of page pointers to gpu
    gpu_memcpy_host2device(gpunlist->plist, gpunlist->plist_h, npages);
    cudaMemset(gpunlist->nbrSizes, 0, particlesPerGpu * sizeof (int));
}


#define XX   0
#define YX   1
#define ZX   2
#define XY   3
#define YY   4
#define ZY   5
#define XZ   6
#define YZ   7
#define ZZ   8
/*
__device__ void  preduceGPU(THREE_MATRIX *h, double *hi, double * x, double * y, double * z)
{
   double da = lround(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
   double db = lround(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));

 *x += h[XX]*da + h[YX]*db;
 *y += h[XY]*da + h[YY]*db;
 *z += h[XZ]*da + h[YZ]*db;
}
 */
/*
__device__ void preduceGPUO3(THREE_MATRIX *h1, double *x, double *y, double *z)
{
   double *h = (double *) h1;
   double hxx,hyy, hzz; 
   hxx = h[XX]; hyy = h[YY]; hzz = h[ZZ]; 

   if (*x > 0.5*hxx ) {  *x += -hxx ; }
   if (*x < -0.5*hxx ) { *x +=  hxx; }

   if (*y > 0.5*hyy ) { *y +=  -hyy ; }
   if (*y < -0.5*hyy ) {*y +=  hyy; }

   if (*z > 0.5*hzz ) { *z +=  -hzz ; }
   if (*z < -0.5*hzz ) {*z +=  hzz; }
}
 */

//__device__ void preduceGPUO3Half(THREE_MATRIX *h1, listtype *x, listtype *y, listtype *z)

__device__ void preduceGPUO3Half(THREE_MATRIX *h1, half *x, half *y, half *z)
{
    double *h = (double *) h1;
    half hxx, hyy, hzz;
    hxx = __float2half(h[XX]);
    hyy = __float2half(h[YY]);
    hzz = __float2half(h[ZZ]);

    if (__hgt(*x, __hmul(0.5, hxx)))
    {
        *x += -hxx;
    }
    if (__hlt(*x, __hmul(-0.5, hxx)))
    {
        *x += hxx;
    }

    if (__hgt(*y, __hmul(0.5, hyy)))
    {
        *y += -hyy;
    }
    if (__hlt(*y, __hmul(-0.5, hyy)))
    {
        *y += hyy;
    }

    if (__hgt(*z, __hmul(0.5, hzz)))
    {
        *z += -hzz;
    }
    if (__hlt(*z, __hmul(-0.5, hzz)))
    {
        *z += hzz;
    }

}

__device__ int modulo(int x, int N)
{
    return (x % N + N) % N;
}

/* Handles converting to half precision */
__global__ void convertToHalf(double *rxbg, double *rybg, double *rzbg, half *rxbg_h, half *rybg_h, half *rzbg_h, int nion)
{
    int pid = blockIdx.x * blockDim.y + threadIdx.y;

    if (pid >= nion)
    {
        return;
    }
    rxbg_h[pid] = __float2half(rxbg[pid]);
    rybg_h[pid] = __float2half(rybg[pid]);
    rzbg_h[pid] = __float2half(rzbg[pid]);
}

__global__ void debugPair(double *rx, double *ry, double *rz, int nion)
{
    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    double rxl = rx[pid];
    double ryl = ry[pid];
    double rzl = rz[pid];
    for (int i = 0; i < nion; i++)
    {
        double dx = rxl - rx[i];
        double dy = ryl - ry[i];
        double dz = rzl - rz[i];
        double d2 = dx * dx + dy * dy + dz*dz;
        if (i != pid && d2 < 9) printf("DebugPair distance i=%d j=%d d2=%f rxl=%f ryl=%f rzl=%f rxg=%f ryg=%f rzg=%f\n", pid, i, d2, rxl, ryl, rzl, rx[i], ry[i], rz[i]);
    }
}

/*
 *Builds neighbor list. rxbg, rybg, rzbf, bin, and plist, and rback  are all bin indexed
 *TODO, it should be possible to exploit symmetry here using a 'diagonal' iteration
 *strategy between 2 bins' interaction matrix, but the atomic needs to be addressed
 */
__global__ void buildList(double *rxbg, double *rybg, double *rzbg,
                          unsigned *bin, int *nbins, int *nbrIdsx, int *nbrIdsy, int *nbrIdsz,
                          int *binHeads, int* numNbrs, int *binCounts, double cut, THREE_MATRIX *h,
                          int **plist, int *nbrSizes, int *rback, int * periodicity, int n, int nion)
{
    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    //int tid = threadIdx.y;
    int tid_x = threadIdx.x;
    //shouldn't have more threads than particles
    if (pid >= nion)
    {
        return;
    }

    //int nbr_count = nbrSizes[pid];

    int ii = rback[pid];
    if (ii >= n) return; //exit it not local
    //double f = 0; 
    //int binsTot = nbins[0]*nbins[1]*nbins[2];
    int local_bin = bin[pid];

    //get local particle assigned to this thread
    double rxl = rxbg[pid];
    double ryl = rybg[pid];
    double rzl = rzbg[pid];

    double cut2 = cut*cut;

    int loc_bin = local_bin; //flattened bin index of current particle pid

    //if bin is on the boundary, some neighboring bins wont exist. skip those bins
    //this shouldn't result in branch diversion if each warp hits same condition
    //use loc_bin to get bin indeces in x,y, and z directions
    int xybins = nbins[0] * nbins[1];
    int loc_bin_z = local_bin / (xybins);
    int loc_bin_xy = (loc_bin - loc_bin_z * (xybins));
    int loc_bin_y = loc_bin_xy / (nbins[0]);
    int loc_bin_x = loc_bin_xy - loc_bin_y * (nbins[0]);

    int *page = plist[0] + PAGE_SIZE*pid;
    //loop over 27 possible neighboring bins
    //numNbrs[0]=3, numNbrs[1]=3, numNbrs[2]=3
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
                //nbrIdsx,nbrIdsy, and nbrIdsz, contain offsets to get bins in contingous x,y, and z directions
                int nbr_bin_z = loc_bin_z + nbrIdsz[k];
                int nbr_bin_y = loc_bin_y + nbrIdsy[j];
                int nbr_bin_x = loc_bin_x + nbrIdsx[i];


                //whether neighbor bin is out of bounds, not considering periodicity
                int out_x = nbr_bin_x < 0 || nbr_bin_x >= nbins[0];
                int out_y = nbr_bin_y < 0 || nbr_bin_y >= nbins[1];
                int out_z = nbr_bin_z < 0 || nbr_bin_z >= nbins[2];

                //whether neighbor bin is out of bounds, but now considering periodicity
                //ie if bin is out of bounds but in periodic direction (ghost bin), it is valid
                //if not periodic, then the neighbor bin is invalid and not considered later
                //if bin is in bounds, it is trivially valid
                int valid_x = (!out_x) || periodicity[0];
                int valid_y = (!out_y) || periodicity[1];
                int valid_z = (!out_z) || periodicity[2];

                //whether bin is a ghost bin
                //a bin is is a ghost bin if it is 'out of bounds' and lies in a periodic direction
                int ghost_x = out_x && periodicity[0];
                int ghost_y = out_y && periodicity[1];
                int ghost_z = out_z && periodicity[2];

                //printf("pid %i ghost_x %i periodic x %i\n", pid, ghost_x, periodicity[0]);
                //skip bins that are not in bounds and not in a periodic direction
                if (!(valid_x && valid_y && valid_z))
                {
                    continue;
                }

                //for debugging, these might go later
                //int nxold = nbr_bin_x;
                //int nyold = nbr_bin_y;
                //int nzold = nbr_bin_z;

                //if bin is a ghost bin, find which bin it is 'replicating'
                if (ghost_x || ghost_y || ghost_z)
                {
                    nbr_bin_x = modulo(nbr_bin_x, nbins[0]);
                    nbr_bin_y = modulo(nbr_bin_y, nbins[1]);
                    nbr_bin_z = modulo(nbr_bin_z, nbins[2]);
                }

                //printf("cur bin %i nbrx %i %i nbry %i %i nbrz %i %i \n",loc_bin, nxold, nbr_bin_x, nyold, nbr_bin_y, nzold, (-1)%3);


                int nbr_bin = nbr_bin_z * nbins[0] * nbins[1] + nbr_bin_y * nbins[0] + nbr_bin_x;
                //now we have to loop over remote particles in neighboring bin.  
                int bin_end = binHeads[nbr_bin] + binCounts[nbr_bin];
                //binHeads points to start of bin. binCounts is length of bin

                int chunk_size = ceil((float) binCounts[nbr_bin] / TPP_NLIST_CONSTRUCT);
                //int chunk_end = min(nbrSizes[pid], (tid_y+1)*chunk_size);

                //thread id flattened across x and y dimensions
                //int tid_flat = blockDim.x*tid+tid_x;
                //int z=blockDim.y*blockDim.x +tid_flat;
                //int block = blockDim.y*blockDim.x;

                for (int jj = binHeads[nbr_bin] + tid_x; jj < binHeads[nbr_bin] + binCounts[nbr_bin]; jj += TPP_NLIST_CONSTRUCT)
                {
                    if (jj >= bin_end)
                    {
                        continue;
                    }
                    if (jj == pid)
                    {
                        continue;
                    } //don't count self interaction

                    double rxj = rxbg[jj];
                    double ryj = rybg[jj];
                    double rzj = rzbg[jj];
                    //printf("x %f y %f z %f\n", rxj, ryj, rzj);
                    //calc partial distances
                    double rxd = rxl - rxj;
                    double ryd = ryl - ryj;
                    double rzd = rzl - rzj;

                    //periodically reduce if processing ghost bin 
                    if (ghost_x || ghost_y || ghost_z)
                    {
                        //double xx= rxd; double yy =ryd; double zz = rzd;
                        PBC(preduceGPUO3(h, &rxd, &ryd, &rzd);)
                            //double *h1 =(double *) h;
                            //printf("x %f - %f = %f ->  %f xx %f yy %f\n",rxl, rxj, xx, rxl-rxd, h1[XX], h1[YY]);
                    }


                    double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
                    //double d= sqrt(d2); //distance
                    //append to neighbor list if distance between particles < cutoff
                    //the tentantive end of the current particle's nlist is nbrSize[pid]. 
                    //thus when we add a particle to the list we need to increment it by 1
                    if (d2 < cut2)
                    {
                        //if(d2<9) printf("BuildList distance pid=%d i=%d j=%d jj=%d d2=%f rxl=%f ryl=%f rzl=%f rxg=%f ryg=%f rzg=%f\n", pid, i, j, jj, d2, rxl, ryl, rzl, rxbg[jj], rybg[jj], rzbg[jj]);

#if(TPP==1)          
                        page[nbr_count] = jj;
                        nbr_count++;
#else
                        int idx = atomicAdd(nbrSizes + pid, 1); //use atomic to append to end of nlist
                        page[idx] = jj; //TODO index into correct page, pages in future not necessarily contiguous   
#endif
                    }
                }
            }
        }
    }
#if(TPP==1)
    nbrSizes[pid] = nbr_count;
#endif

}

__global__ void checkbounds(const int nitems, const int items[], const int maxval, int failcount[1])
{
    const int nt = blockDim.x * blockDim.y;
    const int tid = blockIdx.x * blockDim.y + threadIdx.y;

    for (int i = tid; i < nitems; i += nt)
    {
        if (items[i] > maxval)
            atomicAdd(failcount, 1);
    }
}

__global__ void buildListHalf(half *rxbg, half *rybg, half *rzbg,
                              unsigned *bin, int *nbins, int *nbrIdsx, int *nbrIdsy, int *nbrIdsz,
                              int *binHeads, int* numNbrs, int *binCounts, double cut, THREE_MATRIX *h,
                              int **plist, int *nbrSizes, int *rback, int * periodicity, int n, int nion)
{
    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    //int tid = threadIdx.y;
    int tid_x = threadIdx.x;
    //shouldn't have more threads than particles
    //   if (pid>=nion){return;}   
    if (pid < nion)
    {
        //convert to half
        half hcut = __float2half(cut);
        half hcut2 = __hmul(hcut, hcut);
        int nbr_count = nbrSizes[pid];

        int ii = rback[pid];
        if (ii < n)
        {
            //double f = 0; 
            //int binsTot = nbins[0]*nbins[1]*nbins[2];
            int local_bin = bin[pid];

            //get local particle assigned to this thread
            half rxl = rxbg[pid];
            half ryl = rybg[pid];
            half rzl = rzbg[pid];

            int loc_bin = local_bin;
            //if bin is on the boundary, some neighboring bins wont exist. skip those bins
            //this shouldn't result in branch diversion if each warp hits same condition
            int xybins = nbins[0] * nbins[1];
            int loc_bin_z = local_bin / (xybins);
            int loc_bin_xy = (loc_bin - loc_bin_z * (xybins));
            int loc_bin_y = loc_bin_xy / (nbins[0]);
            int loc_bin_x = loc_bin_xy - loc_bin_y * (nbins[0]);

            int *page = plist[0] + PAGE_SIZE*pid;
            //loop over 27 possible neighboring bins
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


                        //whether neighbor bin is out of bounds, not considering periodicity
                        int out_x = nbr_bin_x < 0 || nbr_bin_x >= nbins[0];
                        int out_y = nbr_bin_y < 0 || nbr_bin_y >= nbins[1];
                        int out_z = nbr_bin_z < 0 || nbr_bin_z >= nbins[2];

                        //whether neighbor bin is out of bounds, but now considering periodicity
                        //ie if bin is out of bounds but in periodic direction (ghost bin), it is valid
                        //if not periodic, then the neighbor bin is invalid and not considered later
                        //if bin is in bounds, it is trivially valid
                        int valid_x = (!out_x) || periodicity[0];
                        int valid_y = (!out_y) || periodicity[1];
                        int valid_z = (!out_z) || periodicity[2];

                        //whether bin is a ghost bin
                        //a bin is is a ghost bin if it is 'out of bounds' and lies in a periodic direction
                        int ghost_x = out_x && periodicity[0];
                        int ghost_y = out_y && periodicity[1];
                        int ghost_z = out_z && periodicity[2];

                        //printf("pid %i ghost_x %i periodic x %i\n", pid, ghost_x, periodicity[0]);
                        //skip bins that are not in bounds and not in a periodic direction
                        if (!(valid_x && valid_y && valid_z))
                        {
                            continue;
                        }

                        //for debugging, these might go later
                        //int nxold = nbr_bin_x;
                        //int nyold = nbr_bin_y;
                        //int nzold = nbr_bin_z;

                        //if bin is a ghost bin, find which bin it is 'replicating'
                        if (ghost_x || ghost_y || ghost_z)
                        {
                            nbr_bin_x = modulo(nbr_bin_x, nbins[0]);
                            nbr_bin_y = modulo(nbr_bin_y, nbins[1]);
                            nbr_bin_z = modulo(nbr_bin_z, nbins[2]);
                        }

                        //printf("cur bin %i nbrx %i %i nbry %i %i nbrz %i %i \n",loc_bin, nxold, nbr_bin_x, nyold, nbr_bin_y, nzold, (-1)%3);


                        int nbr_bin = nbr_bin_z * nbins[0] * nbins[1] + nbr_bin_y * nbins[0] + nbr_bin_x;
                        //now we have to loop over remote particles in neighboring bin.  
                        int bin_end = binHeads[nbr_bin] + binCounts[nbr_bin];
                        //binHeads points to start of bin. binCounts is length of bin

                        int chunk_size = ceil((float) binCounts[nbr_bin] / TPP_NLIST_CONSTRUCT);
                        //int chunk_end = min(nbrSizes[pid], (tid_y+1)*chunk_size);

                        //thread id flattened across x and y dimensions
                        //int tid_flat = blockDim.x*tid+tid_x;
                        //int z=blockDim.y*blockDim.x +tid_flat;
                        //int block = blockDim.y*blockDim.x;

                        for (int jj = binHeads[nbr_bin] + tid_x; jj < binHeads[nbr_bin] + binCounts[nbr_bin]; jj += TPP_NLIST_CONSTRUCT)
                        {
                            if (jj >= bin_end)
                            {
                                continue;
                            }
                            if (jj == pid)
                            {
                                continue;
                            } //don't count self interaction

                            half rxj = rxbg[jj];
                            half ryj = rybg[jj];
                            half rzj = rzbg[jj];
                            //printf("x %f y %f z %f\n", rxj, ryj, rzj);
                            //calc partial distances

                            half rxd = __hsub(rxl, rxj);
                            half ryd = __hsub(ryl, ryj);
                            half rzd = __hsub(rzl, rzj);
                            //periodically reduce if processing ghost bin 
                            if (ghost_x || ghost_y || ghost_z)
                            {
                                //half xx= rxd; half yy =ryd; half zz = rzd;
                                preduceGPUO3Half(h, &rxd, &ryd, &rzd);
                                //double *h1 =(double *) h;
                                //printf("x %f - %f = %f ->  %f xx %f yy %f\n",rxl, rxj, xx, rxl-rxd, h1[XX], h1[YY]);
                            }


                            //distance^2 between particles           
                            half d2 = __hadd(__hmul(rxd, rxd), __hadd(__hmul(ryd, ryd), __hmul(rzd, rzd))); //distance^2 between particles           
                            //append to neighbor list if distance between particles < cutoff
                            //the tentantive end of the current particle's nlist is nbrSize[pid]. 
                            //thus when we add a particle to the list we need to increment it by 1
                            if (__hlt(d2, hcut2))
                            {
                                //				int idx = atomicAdd(nbrSizes+pid, 1); //use atomic to append to end of nlist
                                //				page[idx]=jj;//TODO index into correct page, pages in future not necessarily contiguous
                                page[nbr_count] = jj;
                                nbr_count++;
                            }
                        }
                    }
                }
            }
            nbrSizes[pid] = nbr_count;
        }
    } // end if statements   
}

//fused coulomb - Lennard-jones MARTINI kernel

__device__ thrust::pair<double, double> NKernel2(CharmmLJGPU_PARMS parms, double r, double rinv, double rinv2, double qI, int i, int j, int sij)
{

    //lennard jones piece    
    CharmmLJGPU_PARMS * ljparms = &parms;
    double qJ = ljparms->charge[j];
    double sigma = ljparms->sigma[sij];
    double eps = ljparms->eps[sij];
    double shift = ljparms->shift[sij];
    //double shift =0;
    double ir = rinv;
    double sigma_r = sigma * ir;
    double s2 = sigma_r*sigma_r;
    double s4 = s2*s2;
    double s6 = s4*s2;
    double s12 = s6*s6;
    double valid = r < ljparms->rcut[sij];
    double dvdr = (24.0 * eps * (s6 - 2.0 * s12)) * rinv2;
    double vij = (4.0 * eps * (s12 - s6) + shift);

    //couloumb piece
    //double ke =parms.ke;
    // double krf =parms.krf;
    // double crf = parms.crf;
    // double iepsilon_r=parms.iepsilon_r ;
    double r2 = r*r;
    double epseduo = parms.ke * qI * qJ * parms.iepsilon_r;
    double eij = epseduo * (rinv + parms.krf * r2 - parms.crf);
    double dedr = epseduo * (2 * parms.krf - rinv2 * rinv); // F=-ke*qi*qj/er*[1/r^2 - 2krf*r]r^/r
    thrust::pair<double, double> p((vij + eij) * valid, (dvdr + dedr) * valid);
    //if (r<9) printf("vij=%f eij=%f valid=%f r=%f s6=%f s12=%f shift=%f sigma=%f eps=%f i=%d j=%d\n", vij, eij, valid, r, s6, s12, shift, sigma, eps, i, j);
    //if (i<20 && j<20) printf("vij=%f eij=%f valid=%f r=%f s6=%f s12=%f shift=%f sigma=%f eps=%f i=%d j=%d\n", vij, eij, valid, r, s6, s12, shift, sigma, eps, i, j);
    return p;
}



//Lennard Jones pair.c potential without divisions. The thrust stuff needs to go

__device__ thrust::pair<double, double> NKernel2(LJGPU parms, double r, double rinv, double rinv2, double qI, int i, int j, int sij)
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

//Lennard Jones potential with divisions. The thrust stuff needs to go

__device__ thrust::pair<double, double> NKernel(LJGPU parms, double r, double rinv, double rinv2, int i, int j)
{
    double sigma_r = parms.sig / r;
    double eps = parms.eps;
    double valid = r < parms.rcut;
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


//herein lie all the force evaluation kernels of new and old. 
//thouest shall ascertain that of these contenders, 
//evalList5 reigns supreme in the domain of performance
//and is henceforth the default kernel 

/*
 *Evaluates neighbor list. 
 *r_back stores back pointers from bin indeces to initial indeces
 *fx, fy, and fz are force arrays that do not use the bin indexing
 *e is the energy array, which also is not bin indexed
 *rxbg, rybg, and rzbg are the x,y,z positions of particles,uses bin indeces
 * plist is the list of pointers to each particle's page/s storing nbr indeces (bin indexed)
 * nbrSizes is the number of neighbors for each particle, uses bin indeces
 * parms stores pair parms
 * n is number of particles local to this rank
 */
template <typename K>
__global__ void evalList1(int *rback, double *fx, double *fy, double *fz,
                          double *rxbg, double *rybg, double *rzbg,
                          double *e, int **plist, int *nbrSizes, K * parms, int n, int nion)
{
    __shared__ double shm_e[SHARED_NLIST * TPP_NLIST_EVAL];
    __shared__ double shm_fx[SHARED_NLIST * TPP_NLIST_EVAL];
    __shared__ double shm_fy[SHARED_NLIST * TPP_NLIST_EVAL];
    __shared__ double shm_fz[SHARED_NLIST * TPP_NLIST_EVAL];

    //pid is used to index into the particle "i"
    //each particle "i" has a team of TPP_NLIST threads to process i's neighbors
    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    int tid = threadIdx.y;
    int tid_y = threadIdx.x; //index of thread within particle i's team
    //shm_e[tid_y+tid*TPP_NLIST_EVAL]=0;
    __syncthreads();
    if (pid >= nion)
    {
        return;
    } //exit if out of bonds
    int ii = rback[pid]; //original index of particle i
    if (ii >= n) return; //exit if particle is not local to this process
    int tid_f = tid * TPP_NLIST_EVAL + tid_y; //flattened threadIdx

    //zero the shared memory
    if (tid_f < SHARED_NLIST * TPP_NLIST_EVAL)
    {
        shm_e[tid_f] = 0;
        shm_fx[tid_f] = 0;
        shm_fx[tid_f] = 0;
        shm_fx[tid_f] = 0;
    }
    __syncthreads();

    double ffx = 0;
    double ffy = 0;
    double ffz = 0;
    double ee = 0;
    //get this particle's corresponding page in the neighbor list
    int *page = plist[0] + PAGE_SIZE*pid;
    double rxl = rxbg[pid];
    double ryl = rybg[pid];
    double rzl = rzbg[pid];
    K p = *parms;

    int chunk_size = ceil((float) nbrSizes[pid] / TPP_NLIST_EVAL);
    int chunk_end = min(nbrSizes[pid], (tid_y + 1) * chunk_size);

    //each thread iterates over it's chunk and accumulates forces &energies to a local register
    for (int i = tid_y * chunk_size; i < chunk_end; i++)
    {
        int j = page[i];
        double rxd = rxl - rxbg[j];
        double ryd = ryl - rybg[j];
        double rzd = rzl - rzbg[j];
        double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
        double rinv = rsqrt(d2); //fast inverse distance
        double rinv2 = rinv*rinv; //1/d^2
        double d = rinv*d2; //distance
        int sij = 0;

        thrust::pair<double, double> energy_dvdroverR = NKernel2(p, d, rinv, rinv2, 1.0, pid, j, sij);
        ee += thrust::get<0>(energy_dvdroverR);
        ffx -= thrust::get<1>(energy_dvdroverR) * rxd;
        ffy -= thrust::get<1>(energy_dvdroverR) * ryd;
        ffz -= thrust::get<1>(energy_dvdroverR) * rzd;
    }

    //commit thread's local force&energy accumulator to shmem for reduction
    shm_e[tid * TPP_NLIST_EVAL + tid_y] = ee;
    shm_fx[tid * TPP_NLIST_EVAL + tid_y] = ffx;
    shm_fy[tid * TPP_NLIST_EVAL + tid_y] = ffy;
    shm_fz[tid * TPP_NLIST_EVAL + tid_y] = ffz;
    __syncthreads();
    double eee = 0;
    double fffx = 0;
    double fffy = 0;
    double fffz = 0;
    if (tid_y == 0)
    {
        //TODO: should probably use shuffle sync
#if (TPP_NLIST_EVAL == 4)
        eee += shm_e[tid * TPP_NLIST_EVAL + 0] + shm_e[tid * TPP_NLIST_EVAL + 1] + shm_e[tid * TPP_NLIST_EVAL + 2] + shm_e[tid * TPP_NLIST_EVAL + 3];
        fffx += shm_fx[tid * TPP_NLIST_EVAL + 0] + shm_fx[tid * TPP_NLIST_EVAL + 1] + shm_fx[tid * TPP_NLIST_EVAL + 2] + shm_fx[tid * TPP_NLIST_EVAL + 3];
        fffy += shm_fy[tid * TPP_NLIST_EVAL + 0] + shm_fy[tid * TPP_NLIST_EVAL + 1] + shm_fy[tid * TPP_NLIST_EVAL + 2] + shm_fy[tid * TPP_NLIST_EVAL + 3];
        fffz += shm_fz[tid * TPP_NLIST_EVAL + 0] + shm_fz[tid * TPP_NLIST_EVAL + 1] + shm_fz[tid * TPP_NLIST_EVAL + 2] + shm_fz[tid * TPP_NLIST_EVAL + 3];
#else
        for (int i = 0; i < TPP_NLIST_EVAL; i++)
        {
            eee += shm_e[tid * TPP_NLIST_EVAL + i];
            fffx += shm_fx[tid * TPP_NLIST_EVAL + i];
            fffy += shm_fy[tid * TPP_NLIST_EVAL + i];
            fffz += shm_fz[tid * TPP_NLIST_EVAL + i];
        }
#endif
        //note that we only write to global memory at the end to reduce global mem writes
        e[ii] += .5 * eee; //*.5 because we double count to reduce branch divergence
        fx[ii] += fffx;
        fy[ii] += fffy;
        fz[ii] += fffz;
    }
}

/*
 *Evaluates neighbor list, but stores pages in shared memory
 *for consecutive reads  (not true anymore, now we just
 *modify our thread iteration pattern to get the consecutive reads
 *eg (for int i =tid_y; i<nbrSizes[pid];i+=TPP_NLIST_EVAL)
 *r_back stores back pointers from bin indeces to initial indeces
 *fx, fy, and fz are force arrays that do not use the bin indexing
 *e is the energy array, which also is not bin indexed
 *rxbg, rybg, and rzbg are the x,y,z positions of particles,uses bin indeces
 * plist is the list of pointers to each particle's page/s storing nbr indeces (bin indexed)
 * nbrSizes is the number of neighbors for each particle, uses bin indeces
 * parms stores pair parms
 * n is number of particles local to this rank
also tried using restrict pointers but this didn't help at all

 */
template <typename K>
__global__ void evalList2(int* rback, double* fx, double* fy, double* fz,
                          double* rxbg, double* rybg, double* rzbg, int * sid, double* e,
                          double *sxx, double *syy, double *szz, double *sxy, double*sxz, double *syz,
                          int** plist, int* nbrSizes, K* parms, THREE_MATRIX *h, int n, int nion)
{
    __shared__ double shm_e[SHARED_NLIST * TPP_NLIST_EVAL]; //TODO use synch shuffle for e&f to reserve shmem for page tables
    __shared__ double shm_fx[SHARED_NLIST * TPP_NLIST_EVAL];
    __shared__ double shm_fy[SHARED_NLIST * TPP_NLIST_EVAL];
    __shared__ double shm_fz[SHARED_NLIST * TPP_NLIST_EVAL];
    //__shared__ int shm_page[PAGE_SIZE*SHARED_NLIST]; //todo turn into union

    /*
            int pid = blockIdx.x*blockDim.x + threadIdx.x;
       int tid = threadIdx.x;
            int tid_y = threadIdx.y;
     */
    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    if (pid >= nion)
    {
        return;
    }
    int tid = threadIdx.y;
    int tid_y = threadIdx.x;

    int ii = rback[pid]; //original index of particle i
    if (ii >= n) return; //exit if particle is not local to this process

    //grab some pointers
    //double *sxx = sion->xx;
    //double *syy = sion->yy;
    //double *szz = sion->zz;
    //double *sxy = sion->xy;
    //double *sxz = sion->xz;
    //double *syz = sion->yz;

    shm_e[tid_y + tid * TPP_NLIST_EVAL] = 0;
    __syncthreads();

    int si = sid[pid];
    double ffx = 0;
    double ffy = 0;
    double ffz = 0;
    double ee = 0;
    DVIRIAL(
        double ssxx = 0;
        double ssyy = 0;
        double sszz = 0;
    )
    VIRIAL(
        double ssxy = 0;
        double ssxz = 0;
        double ssyz = 0;
    )

            const int * __restrict__ page = plist[0] + PAGE_SIZE*pid;
    double rxl = rxbg[pid];
    double ryl = rybg[pid];
    double rzl = rzbg[pid];
    double qI = parms->charge[pid];
    K p = *parms;

    //each thread iterates over it's chunk and accumulates forces &energies to a local register
    //printf("ii %i size %i\n", ii, nbrSizes[pid]);
    int chunk_size = ceil((float) nbrSizes[pid] / TPP_NLIST_EVAL);
    int chunk_start = tid_y*chunk_size;
    int chunk_end = min((tid_y + 1) * chunk_size, nbrSizes[pid]);
    for (int i = tid_y; i < nbrSizes[pid]; i += TPP_NLIST_EVAL)
        //for (int i =chunk_start; i<chunk_end;i+=1)
    {
        //int j = shm_page[k*block +tid_flat];
        int j = page[i];
        //k++;
        double rxd = rxbg[j] - rxl;
        double ryd = rybg[j] - ryl;
        double rzd = rzbg[j] - rzl;
        PBC(preduceGPUO3(h, &rxd, &ryd, &rzd);)
                //PBC(preduceGPUO3x(h, rxd,ryd,rzd);)
                double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
        double rinv = rsqrt(d2); //fast inverse  distance 1/sqrt(d^2) = 1/d
        double rinv2 = rinv*rinv; //inverse distance squared 1/d^2
        double d = rinv*d2; //distance
        int sj = sid[j]; //species id


        int sij = sj + parms->nspecies * si;
        thrust::pair<double, double> energy_dvdroverR = NKernel2(p, d, rinv, rinv2, qI, pid, j, sij);
        ee += thrust::get<0>(energy_dvdroverR);
        double fxij = thrust::get<1>(energy_dvdroverR) * rxd;
        double fyij = thrust::get<1>(energy_dvdroverR) * ryd;
        double fzij = thrust::get<1>(energy_dvdroverR) * rzd;

        ffx += fxij;
        ffy += fyij;
        ffz += fzij;
        DVIRIAL(
            ssxx += fxij*rxd;
            ssyy += fyij*ryd;
            sszz += fzij*rzd;
        )
        VIRIAL(
            ssxy += fxij*ryd;
            ssxz += fxij*rzd;
            ssyz += fyij*rzd;
        )
    }


    //start reductions
    //commit thread's local force&energy accumulator to shmem for reduction
    shm_e[tid * TPP_NLIST_EVAL + tid_y] = ee;
    shm_fx[tid * TPP_NLIST_EVAL + tid_y] = ffx;
    shm_fy[tid * TPP_NLIST_EVAL + tid_y] = ffy;
    shm_fz[tid * TPP_NLIST_EVAL + tid_y] = ffz;
    __syncthreads();

    double eee = 0;
    double fffx = 0;
    double fffy = 0;
    double fffz = 0;

    //accumulate each TPP thread's contribution to get total force/energy for a particle
    //unrolling for len 4
    //TODO: use a shuffle sync reduction instead so we can avoid shmem entirely in reduction
    if (tid_y == 0)
    {
        for (int i = 0; i < TPP_NLIST_EVAL; i++)
        {
            eee += shm_e[tid * TPP_NLIST_EVAL + i];
            fffx += shm_fx[tid * TPP_NLIST_EVAL + i];
            fffy += shm_fy[tid * TPP_NLIST_EVAL + i];
            fffz += shm_fz[tid * TPP_NLIST_EVAL + i];
        }
        e[ii] += .5 * eee; //.5 because we double count to reduce branch divergence
        fx[ii] += fffx;
        fy[ii] += fffy;
        fz[ii] += fffz;
    }
    //return;
    __syncthreads();

    DVIRIAL(
            //reduce virials
            double sssxx = 0;
            double sssyy = 0;
            double ssszz = 0;
            //double sssxy=0;
            //double sssxz=0;
            //double sssyz=0;

            shm_fx[tid * TPP_NLIST_EVAL + tid_y] = ssxx;
            shm_fy[tid * TPP_NLIST_EVAL + tid_y] = ssyy;
            shm_fz[tid * TPP_NLIST_EVAL + tid_y] = sszz;
            __syncthreads();
            if (tid_y == 0)
    {
            for (int i = 0; i < TPP_NLIST_EVAL; i++)
        {
            sssxx += shm_fx[tid * TPP_NLIST_EVAL + i];
            sssyy += shm_fy[tid * TPP_NLIST_EVAL + i];
            ssszz += shm_fz[tid * TPP_NLIST_EVAL + i];
        }
            sxx[pid] += .5 * sssxx;
            syy[pid] += .5 * sssyy;
            szz[pid] += .5 * ssszz;
    }
            )
    return;
}

/*
Uses a shuffle sync to do reduction's across a particle's thread team
This indeed uses less shared memory but actually takes hair more time 
more than the shared memory reduction,
 */
template <typename K>
__global__ void evalList3(int* rback, double* fx, double* fy, double* fz,
                          double* rxbg, double* rybg, double* rzbg, int * sid, double* e, double *sxx, double *syy, double *szz, double *sxy, double*sxz, double *syz,
                          int** plist, int* nbrSizes, K* parms, THREE_MATRIX *h, int n, int nion)
{
    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    if (pid >= nion)
    {
        return;
    }
    int tid = threadIdx.y;
    int tid_y = threadIdx.x;
    int ii = rback[pid]; //original index of particle i
    if (ii >= n) return; //exit if particle is not local to this process

    int si = sid[pid];
    double ffx = 0;
    double ffy = 0;
    double ffz = 0;
    double ee = 0;
    DVIRIAL(
            double ssxx = 0;
            double ssyy = 0;
            double sszz = 0;
            )
        VIRIAL(
               double ssxy = 0;
               double ssxz = 0;
               double ssyz = 0;
               )

        const int * __restrict__ page = plist[0] + PAGE_SIZE*pid;
    double rxl = rxbg[pid];
    double ryl = rybg[pid];
    double rzl = rzbg[pid];
    double qI = parms->charge[pid];
    K p = *parms;

    //each thread iterates over it's chunk and accumulates forces &energies to a local register
    for (int i = tid_y; i < nbrSizes[pid]; i += TPP_NLIST_EVAL)
    {
        int j = page[i];
        double rxd = rxbg[j] - rxl;
        double ryd = rybg[j] - ryl;
        double rzd = rzbg[j] - rzl;
        PBC(preduceGPUO3(h, &rxd, &ryd, &rzd);)
            double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
        double rinv = rsqrt(d2); //fast inverse  distance 1/sqrt(d^2) = 1/d
        double rinv2 = rinv*rinv; //inverse distance squared 1/d^2
        double d = rinv*d2; //distance
        int sj = sid[j]; //species id

        int sij = sj + parms->nspecies * si;
        thrust::pair<double, double> energy_dvdroverR = NKernel2(p, d, rinv, rinv2, qI, pid, j, sij);
        ee += thrust::get<0>(energy_dvdroverR);
        double fxij = thrust::get<1>(energy_dvdroverR) * rxd;
        double fyij = thrust::get<1>(energy_dvdroverR) * ryd;
        double fzij = thrust::get<1>(energy_dvdroverR) * rzd;

        ffx += fxij;
        ffy += fyij;
        ffz += fzij;
        DVIRIAL(
                ssxx += fxij*rxd;
                ssyy += fyij*ryd;
                sszz += fzij*rzd;
                )
            VIRIAL(
                   ssxy += fxij*ryd;
                   ssxz += fxij*rzd;
                   ssyz += fyij*rzd;
                   )
    }

    //start reductions
    double eee = ee;
    double fffx = ffx;
    double fffy = ffy;
    double fffz = ffz;
    DVIRIAL(
            double sssxx = ssxx;
            double sssyy = ssyy;
            double ssszz = sszz;
            )
        VIRIAL(
               double sssxy = ssxy;
               double sssxz = ssxz;
               double sssyz = ssyz;
               )

#define MASK 0xff
        //reduce energy, forces, and virials across team of TPP_NLIST_EVAL threads
    for (int i = TPP_NLIST_EVAL / 2; i > 0; i /= 2)
    {
        eee += __shfl_down_sync(MASK, eee, i, TPP_NLIST_EVAL);
        fffx += __shfl_down_sync(MASK, fffx, i, TPP_NLIST_EVAL);
        fffy += __shfl_down_sync(MASK, fffy, i, TPP_NLIST_EVAL);
        fffz += __shfl_down_sync(MASK, fffz, i, TPP_NLIST_EVAL);
        DVIRIAL(
                sssxx += __shfl_down_sync(MASK, sssxx, i, TPP_NLIST_EVAL);
                sssyy += __shfl_down_sync(MASK, sssyy, i, TPP_NLIST_EVAL);
                ssszz += __shfl_down_sync(MASK, ssszz, i, TPP_NLIST_EVAL);
                )
            VIRIAL(
                   sssxy += __shfl_down_sync(MASK, sssxy, i, TPP_NLIST_EVAL);
                   sssxz += __shfl_down_sync(MASK, sssxz, i, TPP_NLIST_EVAL);
                   sssyz += __shfl_down_sync(MASK, sssyz, i, TPP_NLIST_EVAL);
                   )
    }
    /*
       //reduce energy, forces, and virials across team of TPP_NLIST_EVAL threads
       for (int i=TPP_NLIST_EVAL/2; i>0; i/=2)
       {
          eee  += __shfl_down_sync(0xffffffff, eee,  i, TPP_NLIST_EVAL);
          fffx += __shfl_down_sync(0xffffffff, fffx, i, TPP_NLIST_EVAL);
          fffy += __shfl_down_sync(0xffffffff, fffy, i, TPP_NLIST_EVAL);
          fffz += __shfl_down_sync(0xffffffff, fffz, i, TPP_NLIST_EVAL);
          VIRIAL(
             sssxx += __shfl_down_sync(0xffffffff, sssxx,  i, TPP_NLIST_EVAL);
             sssyy += __shfl_down_sync(0xffffffff, sssyy, i, TPP_NLIST_EVAL);
             ssszz += __shfl_down_sync(0xffffffff, ssszz, i, TPP_NLIST_EVAL);
             sssxy += __shfl_down_sync(0xffffffff, sssxy, i, TPP_NLIST_EVAL);
             sssxz += __shfl_down_sync(0xffffffff, sssxz, i, TPP_NLIST_EVAL);
             sssyz += __shfl_down_sync(0xffffffff, sssyz, i, TPP_NLIST_EVAL); 
         )
       }
     */
    //thread zero in the team of TPP_NLIST_EVAL threads has the reduce value, so
    //write that
    if (tid_y == 0)
    {
        //note that we only write to global memory at the end to reduce global mem writes
        e[pid] += .5 * eee; //.5 because we double count to reduce branch divergence
        fx[ii] += fffx;
        fy[ii] += fffy;
        fz[ii] += fffz;
        DVIRIAL(
                sxx[pid] += .5 * sssxx;
                syy[pid] += .5 * sssyy;
                szz[pid] += .5 * ssszz;
                )
            VIRIAL(
                   sxy[pid] += .5 * sssxy;
                   sxz[pid] += .5 * sssxz;
                   syz[pid] += .5 * sssyz;
                   )
    }
}



//in launch bounds, 32 is #threads/block
//14 is the number of blocks we want running simultaneously

template <typename K>
__global__ void
__launch_bounds__(32, 14)
evalList5(const int* rback, double * __restrict__ fx, double* __restrict__ fy, double* __restrict__ fz,
          const double* rxbg, const double* __restrict__ rybg, const double* rzbg, int * sid, double * __restrict__ results, double* e,
          double *sxx, double *syy, double *szz, double *sxy, double*sxz, double *syz,
          int** plist, const int* nbrSizes, const K* parms, const THREE_MATRIX *h, int n, int nion)
{
    //	__shared__ double shm_e[SHARED_NLIST*NUMBUFFS];
    int pid = blockIdx.x * blockDim.y + threadIdx.y;
    if (pid >= nion)
    {
        return;
    }

    //tid_y refers to thread teams/particle
    //int tid = threadIdx.y;
    int tid_y = threadIdx.x;
    int ii = rback[pid]; //original index of particle i
    if (ii >= n) return; //exit if particle is not local to this process

    //the h matrix (the box dimensions)
    double hxx = h->xx;
    double hyy = h->yy;
    double hzz = h->zz;

    //sid is just species id
    int si = sid[pid];
    //DVIRIAL(numResultBuffers=7);

    //temp results holds energy, forces, and then virials
    double tempResults[7] = {0, 0, 0, 0, 0, 0, 0};

    /* 
       double ffx=0;
            double ffy=0;
            double ffz=0;
            double ee=0;
     */
    /*
     DVIRIAL(
        double ssxx=0;
        double ssyy=0;
        double sszz=0;
     )
     VIRIAL(
        double ssxy=0;
        double ssxz=0;
        double ssyz=0;
     )
     */

    //get a pointer to the start of the neighbor list
    //PAGE size is maximum number of neighbors per particle
    const int * __restrict__ page = plist[0] + PAGE_SIZE*pid;

    //rxl, ryl, rzl refere to current particle position
    double rxl = rxbg[pid];
    double ryl = rybg[pid];
    double rzl = rzbg[pid];
    double qI = parms->charge[pid];

    K p = *parms;
    //tells you exactly how many neighbors particle pid has
    int NBRSIZE = nbrSizes[pid];
    //each thread iterates over it's chunk and accumulates forces &energies to a local register
    for (int i = tid_y; i < NBRSIZE; i += TPP_NLIST_EVAL)
    {
        //calculate difference in positions between current particle and neighbor particle
        int j = page[i];
        double rxd = rxbg[j] - rxl;
        double ryd = rybg[j] - ryl;
        double rzd = rzbg[j] - rzl;

        //do a preduce
        PBC(preduceGPUO3x(hxx, hyy, hzz, rxd, ryd, rzd);)

            double d2 = rxd * rxd + ryd * ryd + rzd*rzd; //distance^2 between particles           
        //if(d2<9) printf("EvalList5:distance pid=%d i=%d j=%d d2=%f rxl=%f ryl=%f rzl=%f rxg=%f ryg=%f rzg=%f\n", pid, i, j, d2, rxl, ryl, rzl, rxbg[j], rybg[j], rzbg[j]);
        double rinv = rsqrt(d2); //fast inverse  distance 1/sqrt(d^2) = 1/d
        double rinv2 = rinv*rinv; //inverse distance squared 1/d^2
        double d = rinv*d2; //distance
        int sj = sid[j]; //species id

        //call the force kernel, that returns 2 vals: energy and dvdr_over_r
        int sij = sj + parms->nspecies * si;
        thrust::pair<double, double> energy_dvdroverR = NKernel2(p, d, rinv, rinv2, qI, pid, j, sij);


        //ee+=thrust::get<0>(energy_dvdroverR);
        tempResults[0] += thrust::get<0>(energy_dvdroverR); //energy
        double ffx = thrust::get<1>(energy_dvdroverR) * rxd; //calculate force x direction
        double ffy = thrust::get<1>(energy_dvdroverR) * ryd;
        double ffz = thrust::get<1>(energy_dvdroverR) * rzd;

        //ffx += fxij;
        tempResults[1] += ffx;
        //ffy += fyij;
        tempResults[2] += ffy;
        //ffz += fzij;
        tempResults[3] += ffz;
        DVIRIAL(
                tempResults[4] += ffx*rxd;
                tempResults[5] += ffy*ryd;
                tempResults[6] += ffz*rzd;
                //if(ii==3)printf("j %i rxd %f vir %f \n", rback[j],rxd, fxij*rxd);
                )
            VIRIAL(
                   ssxy += fxij*ryd;
                   ssxz += fxij*rzd;
                   ssyz += fyij*rzd;
                   )
    }

    //start reductions
    /*
         double eee=ee;
         double fffx=ffx;
         double fffy=ffy;
         double fffz=ffz;
     */
    /*   tempResults[0]=.5*ee;
       tempResults[1]=ffx;
       tempResults[2]=ffy;
       tempResults[3]=ffz;
     */ VIRIAL(
                 double sssxx = ssxx;
                 double sssyy = ssyy;
                 double ssszz = sszz;
                 double sssxy = ssxy;
                 double sssxz = ssxz;
                 double sssyz = ssyz;
                 )

        //reduce energy, forces, and virials across team of TPP_NLIST_EVAL threads
    for (int i = TPP_NLIST_EVAL / 2; i > 0; i /= 2)
    {
        /*
        eee  += __shfl_down_sync(0xffffffff, eee,  i, TPP_NLIST_EVAL);
        fffx += __shfl_down_sync(0xffffffff, fffx, i, TPP_NLIST_EVAL);
        fffy += __shfl_down_sync(0xffffffff, fffy, i, TPP_NLIST_EVAL);
        fffz += __shfl_down_sync(0xffffffff, fffz, i, TPP_NLIST_EVAL);
         */
        tempResults[0] += __shfl_down_sync(0xffffffff, tempResults[0], i, TPP_NLIST_EVAL);
        tempResults[1] += __shfl_down_sync(0xffffffff, tempResults[1], i, TPP_NLIST_EVAL);
        tempResults[2] += __shfl_down_sync(0xffffffff, tempResults[2], i, TPP_NLIST_EVAL);
        tempResults[3] += __shfl_down_sync(0xffffffff, tempResults[3], i, TPP_NLIST_EVAL);

        DVIRIAL(
                tempResults[4] += __shfl_down_sync(0xffffffff, tempResults[4], i, TPP_NLIST_EVAL);
                tempResults[5] += __shfl_down_sync(0xffffffff, tempResults[5], i, TPP_NLIST_EVAL);
                tempResults[6] += __shfl_down_sync(0xffffffff, tempResults[6], i, TPP_NLIST_EVAL);
                )
            VIRIAL(
                   sssxy += __shfl_down_sync(0xffffffff, sssxy, i, TPP_NLIST_EVAL);
                   sssxz += __shfl_down_sync(0xffffffff, sssxz, i, TPP_NLIST_EVAL);
                   sssyz += __shfl_down_sync(0xffffffff, sssyz, i, TPP_NLIST_EVAL);
                   )
    }


    if (tid_y == 0)
    {
        int start = NUMBUFFS*ii;
        //note that we only write to global memory at the end to reduce global mem writes
        results[start + 0] += .5 * tempResults[0]; //.5 because we double count to reduce branch divergence
        results[start + 1] += tempResults[1];
        results[start + 2] += tempResults[2];
        results[start + 3] += tempResults[3];
        DVIRIAL(
                results[start + 4] += tempResults[4];
                results[start + 5] += tempResults[5];
                results[start + 6] += tempResults[6];
        )
        VIRIAL(
               sxy[pid] += .5 * sssxy;
               sxz[pid] += .5 * sssxz;
               syz[pid] += .5 * sssyz;
        )
    }

    /*
       if (tid_y==0)
       {
         shm_e[tid*NUMBUFFS+ 0]=.5*tempResults[0];
         shm_e[tid*NUMBUFFS+1]=tempResults[1];
         shm_e[tid*NUMBUFFS+2]=tempResults[2];
         shm_e[tid*NUMBUFFS+3]=tempResults[3];
         DVIRIAL(
            shm_e[tid*NUMBUFFS+4]=.5*tempResults[4];
            shm_e[tid*NUMBUFFS+5]=.5*tempResults[5];
            shm_e[tid*NUMBUFFS+6]=.5*tempResults[6];
         )
          //if (ii<5){printf("gpu m %i %f\n", ii, tempResults[4]);}
         //shm_e[tid*5+tid_y]=tempResults[tid_y];
       }
     */
    /*
       __syncthreads();
       if (tid_y<NUMBUFFS)
       {
                    int start = NUMBUFFS*ii;
                    results[start+tid_y]+=shm_e[tid*NUMBUFFS+tid_y]; //.5 because we double count to reduce branch divergence
       }
     */
    /*
       if (tid_y==0)
       {
       int numResultBuffers=4;
       VIRIAL(numResultBuffers=10;)
       int start = numResultBuffers*pid;
            //note that we only write to global memory at the end to reduce global mem writes
            results[start]+=.5*eee; //.5 because we double count to reduce branch divergence
            results[start+1]+=fffx;
            results[start+2]+=fffy;
            results[start+3]+=fffz;
       VIRIAL(
          results[start+4]+=.5*sssxx;
          results[start+5]+=.5*sssyy;
          results[start+6]+=.5*ssszz;
          results[start+7]+=.5*sssxy;
          results[start+8]+=.5*sssxz;
          results[start+9]+=.5*sssyz;
       )
       }
     */
}

/*
{
        //note that we only write to global memory at the end to reduce global mem writes
        results[pid]+=.5*eee; //.5 because we double count to reduce branch divergence
        results[nion+pid]+=fffx;
        results[2*nion+pid]+=fffy;
        results[3*nion+pid]+=fffz;
   VIRIAL(
      results[4*nion+pid]+=.5*sssxx;
      results[5*nion+pid]+=.5*sssyy;
      results[6*nion+pid]+=.5*ssszz;
      results[7*nion+pid]+=.5*sssxy;
      results[8*nion+pid]+=.5*sssxz;
      results[9*nion+pid]+=.5*sssyz;
   )
   }




 */


__global__ void unpermuteForEvalList5(int* rback, double* fx, double* fy, double* fz,
                                      double *results, double* e, double *sxx, double *syy, double *szz,
                                      double *sxy, double*sxz, double *syz, int n, int nion)
{
    //int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int ii = pid; //original index of particle i
    if (pid >= nion)
    {
        return;
    }
    if (pid >= n) return; //exit if particle is not local to this process

    int start = NUMBUFFS*pid; //*rback[pid];
    e[ii] = results[start];

    fx[ii] = results[start + 1];
    fy[ii] = results[start + 2];
    fz[ii] = results[start + 3];


    DVIRIAL(
            //sxx[ii] = results[start+4];
            //syy[ii] = results[start+5];
            //szz[ii] = results[start+6];
            atomicAdd(sxx + ii, results[start + 4]);
            atomicAdd(syy + ii, results[start + 5]);
            atomicAdd(szz + ii, results[start + 6]);
            //  sxy[ii] = results[start+7];
            //  sxz[ii] = results[start+8];
            //  syz[ii] = results[start+9];
            )

    return;
}

template <typename T>
__global__ void unpermuteEnergyForEvalList5(int* rback, T energies, T sxxs, T syys, T szzs, T results, int n, int nion)
{
    //int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int ii = pid; //original index of particle i
    if (pid >= nion)
    {
        return;
    }
    if (pid >= n) return; //exit if particle is not local to this process

    int start = NUMBUFFS*pid; //*rback[pid];
    energies[ii] = results[start];

    DVIRIAL(
            //sxxs[ii] =  results[start+4];
            //syys[ii] = results[start+5];
            //szzs[ii] = results[start+6];
            atomicAdd(sxxs + ii, .5 * results[start + 4]);
            atomicAdd(syys + ii, .5 * results[start + 5]);
            atomicAdd(szzs + ii, .5 * results[start + 6]);
            )

    return;
}

/*
Utility function for printing a neighbor table. 
 */
void printNTable(SYSTEM *sys) 
{
    ACCELERATOR * accelerator = accelerator_getAccelerator(NULL);
    GPUCUDAPARMS *accParms = (GPUCUDAPARMS*) accelerator->parms; 
    //grab some values/pointers for convenience
    GPUNLIST *gnlist=sys->collection->gnlist;

    int n = sys->collection->state->nlocal;
    //STATE *gsh = sys->collection->gpustate_h;
    int nion = sys->collection->state->nion;

    //allocate space for nlist on cpu
    int maxPagesPerParticle=accParms->maxPagesPerParticle;
    
    int * mem = (int *) malloc(sizeof (int)*maxPagesPerParticle * PAGE_SIZE * nnppgpu);
    int *nbrsize = (int *) malloc(sizeof (int)*n);
    int *rback = (int*) malloc(sizeof (int)*nion);

    //memcopy nlist and other relevant metadata to cpuy
    gpu_memcpy_device2host(rback, gnlist->r_backg, nion);
    gpu_memcpy_device2host(nbrsize, gnlist->nbrSizes, n);
    gpu_memcpy_device2host(mem, gnlist->mem, maxPagesPerParticle * PAGE_SIZE * nnppgpu);


    for (int i = 0; i < n; i++)
        printf("r %i %i\n", i, rback[i]);
    for (int i = 0; i < n; i++)
        printf("len %i %i\n", rback[i], nbrsize[rback[i]]);

    //print list
    for (int i = 0; i < n; i++)
    {
        printf("  particle %i\n  ", rback[i]);
        for (int j = 0; j < nbrsize[i]; j++)
        {
            printf(" %i ", rback[mem[PAGE_SIZE * i + j]]);
        }
        printf("\n");
    }
}

/*Constructs Neigbor List on GPU
 *note that 'rcut' is actually the pair cutoff+skin layer
 *where the skin is the extra "buffer" space added to
cutoff radius that allow neighbor list to stay valid
for multiple time steps 
a*/
void constructList(SYSTEM *sys, double rcut) 
{
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    ACCELERATOR * accelerator = accelerator_getAccelerator(NULL);
    GPUCUDAPARMS *accParms = (GPUCUDAPARMS*) accelerator->parms; 
    //grab some pointers for convenience
    GPUNLIST *gnlist=sys->collection->gnlist;
    int n = sys->collection->state->nlocal;
    int nion = sys->collection->state->nion;
    STATE *gsh = sys->collection->gpustate_h;
    // printf("rebuilding list \n");
    //printf("nlocal %i nion %i\n", n, nion);
    //size_t free, total;
    //CUDA_SAFE_CALL(cudaMemGetInfo	(&free ,	&total));
    //printf("used %zu free %zu\n", (total-free)/1048576, free/1048576);


    CUDA_SAFE_CALL(cudaPeekAtLastError());
    cudaMemset(gnlist->rxbg, 0, sizeof (double)*nion);
    cudaMemset(gnlist->rybg, 0, sizeof (double)*nion);
    cudaMemset(gnlist->rzbg, 0, sizeof (double)*nion);
    cudaMemset(gsh->fx, 0, sizeof (double)*nion);
    cudaMemset(gsh->fy, 0, sizeof (double)*nion);
    cudaMemset(gsh->fz, 0, sizeof (double)*nion);
    cudaMemset(gnlist->nbrSizes, 0, sizeof (int)*nion);

    CUDA_SAFE_CALL(cudaPeekAtLastError());
    //memcopy state from cpu to gpu
    SIMULATE *sim = (SIMULATE*) (sys->parent);
    //DDC* ddc=sim->ddc;
    int gpu_integrate = sim->integrator->uses_gpu;
    if (gpu_integrate == 0)
    {
        sendGPUState(sys, sys->nion);
    }

    //sort particles into bins
    binParticlesGPU(sys, rcut);
    /*
        int bsize=256;
        int gsize=ceil((float) nion / bsize);
        printf("Debug Pair orig\n\n");

        debugPair<<<gsize, bsize>>>(gsh->rx, gsh->ry, gsh->rz, nion);
    
        printf("Debug Pair new\n\n");
        debugPair<<<gsize, bsize>>>(gnlist->rxbg, gnlist->rybg, gnlist->rzbg, nion);
     */
    //launch neighbor list force/energy kernel
    int blockSize = 64;
    int gridSize = ceil((float) nion / blockSize);
    dim3 blockSize3(TPP_NLIST_CONSTRUCT, blockSize);

    //launch neighbor list construction kernel
    buildList << <gridSize, blockSize3>>>(gnlist->rxbg, gnlist->rybg, gnlist->rzbg,
        gnlist->binsbg, gnlist->nbinsg, gnlist->nbrIdsx, gnlist->nbrIdsy, gnlist->nbrIdsz,
        gnlist->binHeadsg, gnlist->numNbrsxyz, gnlist->binCountsg2, rcut, gnlist->hmat_g,
        gnlist->plist, gnlist->nbrSizes, gnlist->r_backg, gnlist->pbc_g, n, nion);

    //checkbonds(const int nitems,const int items[],const int maxval,int failcount[1])
    if(accParms->checkBounds == 1)
    {
        int maxPagesPerParticle=accParms->maxPagesPerParticle;
        
        int failcount = 0, *failcount_ptr = &failcount, *g_failcount;
        gpu_allocator(g_failcount, 1);
        gpu_memcpy_host2device(g_failcount, failcount_ptr, 1);
        checkbounds << <gridSize, blockSize3>>>(nion, gnlist->nbrSizes, PAGE_SIZE * maxPagesPerParticle, g_failcount);
        gpu_memcpy_device2host(failcount_ptr, g_failcount, 1);
        assert(failcount == 0 || "Error, not enough neighboe spapce per particle" == NULL);
        cudaFree(g_failcount);
    }

    /*
            //launch neighbor list construction kernel in half precision
       convertToHalf<<<gridSize, blockSize3>>>(gnlist->rxbg, gnlist->rybg, gnlist->rzbg, gnlist->gpu_types->rxbg_h, gnlist->gpu_types->rybg_h, gnlist->gpu_types->rzbg_h, nion);
            buildListHalf<<<gridSize, blockSize3>>>(gnlist->gpu_types->rxbg_h, gnlist->gpu_types->rybg_h, gnlist->gpu_types->rzbg_h,
                              gnlist->binsbg, gnlist->nbinsg, gnlist->nbrIdsx, gnlist->nbrIdsy, gnlist->nbrIdsz, 
                              gnlist->binHeadsg, gnlist->numNbrsxyz, gnlist->binCountsg2,  rcut,gnlist->hmat_g,
                            gnl->plist, gnl->nbrSizes,gnlist->r_backg,gnlist->pbc_g, n, nion);
     */
    CUDA_SAFE_CALL(cudaPeekAtLastError());
    //printNTable(sys);	exit(0); //print neighbor table

}

/*Launches pair/force evaluation kernel on GPU
 *function also memcopies pair info to gpu...this
 *needs to be fixed. Note that 'e' currently
 *is unaffected by this func, and is only changed
 *later during the sendForceEnergyToHost function call
 *where collection->e_all gets reduced and stored in 'e'
 *this may change, hence reason for leaving the arg here
 */
template <typename K>
void evalList(SYSTEM *sys, K* parms, int *sid, ETYPE *e) 
{
    //grab some values/pointers for convenience
    SIMULATE *sim = (SIMULATE*) (sys->parent);
    GPUNLIST *gnlist=sys->collection->gnlist;
    COLLECTION *collection = sys->collection;
    int n = collection->state->nlocal;
    int nion = collection->state->nion;
    STATE *gsh = collection->gpustate_h;

    //launch neighbor list force/energy kernel
    int blockSize = SHARED_NLIST;

    int gridSize = ceil((float) nion / blockSize);
    dim3 blockSize3(TPP_NLIST_EVAL, blockSize);
    //printf("evalList5 gridSize=%d  bock %d  %d\n", gridSize, TPP_NLIST_EVAL, blockSize); 
    evalList5<K> << <gridSize, blockSize3>>>(gnlist->r_backg, gsh->fx, gsh->fy, gsh->fz,
        gnlist->rxbg, gnlist->rybg, gnlist->rzbg, sid, gnlist->results, gnlist->e_all,
        gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy,
        gnlist->gpu_sion->xz, gnlist->gpu_sion->yz, gnlist->plist, gnlist->nbrSizes, parms, gnlist->hmat_g, n, nion);

    //int DSIZE=collection->state->nion*7*sizeof(double);
    //double *h_data = (double *)malloc(DSIZE);
    //cudaMemcpy(h_data, gnl->results, DSIZE, cudaMemcpyDeviceToHost);
    //for(int i=0; i<collection->state->nion; i++){
    //if(h_data[i*7]>0.1) printf("i=%d  evalList5 e_all=%f\n", i, h_data[i*7]);
    //}


    unpermuteEnergyForEvalList5 << <gridSize, blockSize>>>(gnlist->listIdsg, gnlist->e_all, gnlist->gpu_sion->xx,
        gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->results, n, nion);

    //DSIZE=collection->state->nion*sizeof(double);
    //cudaMemcpy(h_data, collection->e_all, DSIZE, cudaMemcpyDeviceToHost);
    //for(int i=0; i<collection->state->nion; i++){
    //if(h_data[i]>0.1) printf("i=%d  unpermuteEnergy e_all=%f\n", i, h_data[i]);
    // }
    /*
   __global__ void evalList2(int*  rback, double* fx, double* fy, double*  fz, 
                             double* rxbg, double*  rybg, double*  rzbg, int * sid, double* e, 
                             double *sxx, double *syy, double *szz, double *sxy, double*sxz, double *syz,
                             int**  plist,  int* nbrSizes, K* parms, THREE_MATRIX *h, int n, int nion)
     */
    /*
            evalList3<K><<<gridSize, blockSize3>>>(gnlist->r_backg, gsh->fx,gsh->fy, gsh->fz,
                                                       gnlist->rxbg, gnlist->rybg, gnlist->rzbg, sid, gnlist->e_all, 
                         gnlist->gpu_sion->xx, gnlist->gpu_sion->yy, gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, 
                         gnlist->gpu_sion->xz, gnlist->gpu_sion->yz,  gnlist->plist, gnlist->nbrSizes, parms,gnlist->hmat_g, n,nion);
     */

    //DDC* ddc=sim->ddc;

    int gpu_integrate = sim->integrator->uses_gpu;

    //if using cpu to integrate
    //put gpu force data in cpu data format for cpu integration
    if (!gpu_integrate)
    {
        unpermuteForEvalList5 << <gridSize, blockSize>>>(gnlist->listIdsg, gsh->fx, gsh->fy, gsh->fz,
            gnlist->results, gnlist->e_all, gnlist->gpu_sion->xx, gnlist->gpu_sion->yy,
            gnlist->gpu_sion->zz, gnlist->gpu_sion->xy, gnlist->gpu_sion->xz, gnlist->gpu_sion->yz,
            n, nion);
    }
    CUDA_SAFE_CALL(cudaPeekAtLastError());
}

/*
This kernel uses shuffle syn reduction instead of shared memory reductions
Therefore, the thread and block layouts are transposed; a particle's team
in the force evaluation now exists in the x thread dimension, and particles
span the y dimension
 */
template <typename K>
void evalList_shuffle(SYSTEM *sys, K* parms, ETYPE *e) 
{
    //grab some values/pointers for convenience
    GPUNLIST *gnlist=sys->collection->gnlist;
    int n = sys->collection->state->nlocal;
    int nion = sys->collection->state->nion;
    STATE *gsh = sys->collection->gpustate_h;

    //launch neighbor list force/energy kernel
    int blockSize = SHARED_NLIST;
    dim3 gridSize = (ceil((float) nion / blockSize));
    dim3 blockSize3(TPP_NLIST_EVAL, blockSize);
    //printf("blocks %i \n", gridSize);
    //printf("shared %i \n", PAGE_SIZE*SHARED_NLIST*32+SHARED_NLIST*TPP_NLIST_EVAL*64);
    evalList1<K> << <gridSize, blockSize3>>>(gnlist->r_backg, gsh->fx, gsh->fy, gsh->fz,
        gnlist->rxbg, gnlist->rybg, gnlist->rzbg,
        gnlist->e_all, gnlist->plist, gnlist->nbrSizes, parms, n, nion);
    CUDA_SAFE_CALL(cudaPeekAtLastError());
}


//used to zero force buffers, energy buffers, velocity buffers  a results buffer (which contains the forces and
//energies and virials, which then get transferred to aforementioned force and energy buffers)
//we initially used memsets, but they were slow because you would call memset many times

__global__ void zeroKernel(STATE * gpu_state, double *e1, double *e_all, double *results,
        double *vx, double *vy, double *vz, int nlocal) 
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= nlocal)return;
    gpu_state->fx[pid] = 0;
    gpu_state->fy[pid] = 0;
    gpu_state->fz[pid] = 0;
    e1[pid] = 0;
    e_all[pid] = 0;
    results[7 * pid] = 0;
    results[7 * pid + 1] = 0;
    results[7 * pid + 2] = 0;
    results[7 * pid + 3] = 0;
    results[7 * pid + 4] = 0;
    results[7 * pid + 5] = 0;
    results[7 * pid + 6] = 0;
    vx[pid] = 0;
    vy[pid] = 0;
    vz[pid] = 0;
}

void zeroGPUForceEnergyBuffers(SYSTEM *sys) 
{
    GPUNLIST *gnlist=sys->collection->gnlist;
    //STATE *gsh = sys->collection->gpustate_h;
    GPUVIRIALS *virials = gnlist->gpu_sion;
    //zero the energy and force  arrays on gpu
    size_t blockSize = 32;
    size_t gridSize = (blockSize + sys->nlocal - 1) / blockSize;
    zeroKernel<<<gridSize, blockSize>>>(sys->collection->gpustate, gnlist->e1, gnlist->e_all, gnlist->results,
            virials->xx, virials->yy, virials->zz, sys->nlocal);
}

void zeroGPUForceEnergyBuffers2(SYSTEM *sys) 
{
    GPUNLIST *gnlist=sys->collection->gnlist;
    //STATE *gsh = sys->collection->gpustate_h;
    //zero the energy and force  arrays on gpu
    cudaMemset(gnlist->e_all, 0, sizeof (double)*sys->nlocal);
    cudaMemset(gnlist->e1, 0, sizeof (double)*sys->nlocal);
    cudaMemset(sys->collection->gpustate_h->fx, 0, sizeof (double)*sys->nlocal);
    cudaMemset(sys->collection->gpustate_h->fy, 0, sizeof (double)*sys->nlocal);
    cudaMemset(sys->collection->gpustate_h->fz, 0, sizeof (double)*sys->nlocal);

    int numResultBuffers = 7;
    //VIRIAL(numResultBuffers=10);
    cudaMemset(gnlist->results, 0, numResultBuffers * sizeof (double)*sys->nion);


    //GPUVIRIALS *virials = sys->collection->gpu_sion_h;
    GPUVIRIALS *virials = gnlist->gpu_sion;
    cudaMemset(virials->xx, 0, sizeof (double)*sys->nlocal);
    cudaMemset(virials->yy, 0, sizeof (double)*sys->nlocal);
    cudaMemset(virials->zz, 0, sizeof (double)*sys->nlocal);
    cudaMemset(virials->xy, 0, sizeof (double)*sys->nlocal);
    cudaMemset(virials->xz, 0, sizeof (double)*sys->nlocal);
    cudaMemset(virials->yz, 0, sizeof (double)*sys->nlocal);
}

//wrapper func for doing pair processing with pair.c LJ potential

void pairProcessNListGpu(SYSTEM *sys, PAIR_PARMS *parms, ETYPE *e) 
{
    GPUNLIST *gnlist=sys->collection->gnlist;
    //STATE * gsh = sys->collection->gpustate_h;
    //printf("processing\n");
    //send particles from host to gpu
    SIMULATE *sim = (SIMULATE*) (sys->parent);
    int gpu_integrate = sim->integrator->uses_gpu;
    if (gpu_integrate == 0)sendGPUState(sys, sys->nion);

    //set force and energy buffers on gpu to zero
    //TODO: this is obviously a huge red flag if you want to 
    //use other potential in addition. you'll have to comment this out
    // in that case
    //zeroGPUForceEnergyBuffers(sys);

    //apply saved permutation (computed in most recent binsort) to quickly re-bin particles received from host
    //permuteParticles(gnlist->listIdsg, gsh->rx, gsh->ry, gsh->rz,gnlist->r_backg,
    //                 gnlist->rxbg, gnlist->rybg, gnlist->rzbg, gnlist->binsg, gnlist->binsbg, sys->nion);
    permuteParticles2(sys->collection, sys->nion);
    //send pair/potential info to gpu
    LJGPU *lj = (LJGPU*) malloc(sizeof (LJGPU));
    LJ_PARMS** ljp = (LJ_PARMS**) parms->parms;
    lj->eps = ljp[0]->eps;
    lj->sig = ljp[0]->sigma;
    lj->shift = ljp[0]->shift;
    lj->rcut = parms->rcut[0].value;
    lj->nspecies = sys->nspecies;
    lj->charge = gnlist->charge_bg;
    LJGPU *ljgpu;
    gpu_allocator(ljgpu, 1);
    gpu_memcpy_host2device(ljgpu, lj, 1);

    //evaluate forces/energy on gpu using gpu neighbor list
    //printf("evaling list \n");
    evalList(sys, ljgpu, gnlist->species_bg, e);
    //send energy and forces from gpu to hsot
    if (gpu_integrate == 0)sendForceEnergyToHost(sys, e);
}


//wrapper func for doing pair processing with charmm/martini potential

void charmmPairGPU(SYSTEM *sys, CHARMMPOT_PARMS *parms, ETYPE *e) 
{
    COLLECTION * collection = sys->collection;
    //STATE * gsh = sys->collection->gpustate_h;

    //copy species data
    STATE *state = sys->collection->state;

    CHARMM_PARMS *charmmParms = parms->charmmParms;
    //int nspecies = parms->nspecies; // has been set to mmff->nAtomType

    //zeroGPUForceEnergyBuffers(sys);
    SIMULATE* simulation = (SIMULATE*) sys->parent;
    DDC* ddc = simulation->ddc;
    PUSH_RANGE("gen species map", 0);

    if (ddc->lastUpdate == sys->loop)
    { // re-assign map upon domain changes

        int speciesIndex[sys->nspecies];
        int *speciesIndexGPU;
        //todo: put this malloc in gpuMemUtils.cu so we don't get memory leak when doing redomains
        gpu_allocator(speciesIndexGPU, sys->nspecies);
        for (int j = 0; j < sys->nspecies; j++)
        {
            speciesIndex[sys->species[j]->index] = getCGLJindexbySpecie(sys->species[j], charmmParms);
        }

        //int sIndex[sys->nion];
        int *sIndex = (int *) ddcMalloc(sys->nion * sizeof (int));
        //unsigned sIndex_blk;
        //int* sIndex = (int *)heapGet(&sIndex_blk);
        //heapEndBlock(sIndex_blk, sys->nion*sizeof(int));
        //needs to get permute later
        for (unsigned i = 0; i < sys->nion; i++)
            sIndex[i] = speciesIndex[state->species[i]->index];

        //int *sIndexPtr = sIndex;
        gpu_memcpy_host2device(parms->gpu_ljparms_h->cg_species_index, sIndex, sys->nion);

        //heapFree(sIndex_blk);
        ddcFree(sIndex);
    }
    POP_RANGE();

    //apply saved permutation (computed in most recent binsort) to quickly re-bin particles received from host
    permuteParticlesCharmm(collection, parms, sys->nion);
    CUDA_SAFE_CALL(cudaPeekAtLastError();)

    //printNTable(sys);  
    //do force calculation
    evalList(sys, parms->gpu_ljparms, parms->gpu_ljparms_h->cg_species_index_b, e);
    PUSH_RANGE("self coulomb sum", 1);
    double q2 = 0.0;
    for (unsigned i = 0; i < sys->nlocal; i++)
    {
        q2 += SQ(state->q[i]); // Should be in the sys->local or sys->nion?
    }

    double keR = ke / parms->epsilon_r;
    POP_RANGE();

    double vEle = -0.5 * q2 * keR * parms->crf; //Self electronic energy   
    e->eion += vEle;
}

