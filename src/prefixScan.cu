#include "prefixScan.h"
#include <cuda_runtime.h>
#include "cudaUtils.h" 
#include <limits.h>
#include <float.h>
#include <sys/time.h>

#define CHECK(x) 


typedef double(*Op)(double, double);

//may neeed to use __forceinline__  here

template <typename T> __device__ T add(T a, T b) {
    return a + b;
}
//naive prefix scan

template <Op op, typename T> __global__ void scan1(T *g_idata, T *g_odata, T *sums, int n) {
    // each thread loads one element from global to shared mem
    int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;

    //load stuff into shared mem.
    //we allocate 2wice as much shmem to create 'negative' buffer space
    __shared__ T sdata[SHARED_PREFIX_SCAN];
    sdata[tid] = 0;
    __syncthreads();
    if (pid < n)
        sdata[tid] = g_idata[pid];
    __syncthreads();
    // do scan in shared mem
    int j = 1;
    T x = 0;
    while (j < blockDim.x) {
        if (tid >= j)
            x = op(sdata[tid], sdata[tid - j]);
        else
            x = sdata[tid];
        __syncthreads();
        sdata[tid] = x;
        __syncthreads();
        j *= 2;
    }
    //commit results to global mem
    if (pid < n) {
        if (tid > 0) {
            g_odata[pid] = sdata[tid - 1];
        } else {
            sums[blockIdx.x] = sdata[blockDim.x - 1];
            g_odata[pid] = 0;
        }
    }
}

//slightly better prefix scan
//previously we used if statement to decide if thread tid needs to
//add arr[tid] and arr[tid-j] (because if j>tid, then no add is needed/possible since tid-j<0)
//here we always do the add by allocating double the shared memory needed,d
//putting the data of interest in the latter half of the shared mem
//and setting the first 1/2 to zero. this prevents branch divergence by removing if statement

template <Op op, typename T> __global__ void scan2(T *g_idata, T *g_odata, T *sums, int n) {
    // each thread loads one element from global to shared mem
    int tid = threadIdx.x;
    int pid = blockIdx.x * blockDim.x + threadIdx.x;

    //load stuff into shared mem.
    //we allocate 2wice as much shmem to create 'negative' buffer space
    __shared__ T sdata[2 * SHARED_PREFIX_SCAN];
    sdata[tid] = 0;
    sdata[tid + SHARED_PREFIX_SCAN] = 0;
    __syncthreads();

    if (pid < n)
        sdata[SHARED_PREFIX_SCAN + tid] = g_idata[pid];
    __syncthreads();

    // do reduction in shared mem
    int j = 1;
    T x = 0;
    while (j < blockDim.x) {
        x = op(sdata[SHARED_PREFIX_SCAN + tid], sdata[SHARED_PREFIX_SCAN + tid - j]);
        __syncthreads();
        sdata[SHARED_PREFIX_SCAN + tid] = x;
        __syncthreads();
        j *= 2;
    }

    //commit results to global mem
    if (pid < n) {
        if (tid > 0)g_odata[pid] = sdata[SHARED_PREFIX_SCAN + tid - 1];
        else {
            sums[blockIdx.x] = sdata[SHARED_PREFIX_SCAN + blockDim.x - 1];
            g_odata[pid] = 0;
        }
    }
}

/*
//like scan2, but avoids 1 of the syncthreads in the while loop
//allocate 4x the shared memory
//than used in scan 1. We allocate 2x the shared mem using the same strategy as in scan 2, and also another 2* the shared mem so we can double buffer when loading and then
//storing in loop, so we don't have to synch between load and store
template <Op op> __global__ void scan3(double *g_idata, double *g_odata,double *sums, int n)
{
   // each thread loads one element from global to shared mem
   int tid = threadIdx.x;
   int pid = blockIdx.x*blockDim.x + threadIdx.x;

   //load stuff into shared mem.
   //we allocate 4x as much shmem to create 'negative' buffer space
   //and to double buffer
   
   sdata[tid]=0;   sdata[tid+3*SHARED_PREFIX_SCAN]=0;
   sdata[tid+2*SHARED_PREFIX_SCAN]=0;
   sdata[tid+1*SHARED_PREFIX_SCAN]=0;
   __syncthreads();
   int b_in =1;
   int b_out 1-b_in;
   int in_start = b_in*2*SHARED_PREFIX_SCAN+SHARED_PREFIX_SCAN;
   int out_start = b_out*2*SHARED_PREFIX_SCAN+SHARED_PREFIX_SCAN;
   if( pid<n)
      sdata[in_start+tid]=g_idata[pid];
   __syncthreads();

   // do reduction in shared mem
   int j = 1;
   double x=0; 
   while(j<blockDim.x) {
      sdata[out_start+tid]=op(sdata[in_start+tid], sdata[in_start+tid - j]);
      __syncthreads();
      j*=2;
      b_in =1-b_in;
      b_out = 1-b_out;
      n_start = b_in*2*SHARED_PREFIX_SCAN+SHARED_PREFIX_SCAN;
      out_start = b_out*2*SHARED_PREFIX_SCAN+SHARED_PREFIX_SCAN;
   }
   
   //commit results to global mem
   if (pid<n)
   {
      if(tid>0 )g_odata[pid] = sdata[SHARED_PREFIX_SCAN+tid-1];
      else{
         sums[blockIdx.x]=sdata[SHARED_PREFIX_SCAN+blockDim.x-1];
         g_odata[pid] = 0;
      }
   }
}
 */

//applies the inter block sum prefix scan to the intra block prefix scan

template <Op op, typename T> __global__ void applySum(T *g_idata, T *g_odata, T *sums, int n) {
    //int tid = threadIdx.x;
    int pid = blockDim.x * blockIdx.x + threadIdx.x;
    if (pid >= n) return;
    T sum = sums[blockIdx.x];
    g_odata[pid] = g_idata[pid] + sum;

}

template <typename T> void launchPlusScan(T *in, T *out, T *sums, T *sums2, int n) {
    int blockSize = min(n, SHARED_PREFIX_SCAN);
    int gridSize = ceil((float) n / blockSize);
    scan1<add><<<gridSize, blockSize>>>(in, out, sums, n);
    scan1<add><<<1, gridSize>>>(sums, sums2, in, gridSize);
    applySum<add><<<gridSize, blockSize>>>(out, out, sums2, n);
}

//note this works for arrays with up to blocksize^2 elems

void plusScan(double *in, double *out, double *sums, double *sums2, int n) {
    int blockSize = min(n, SHARED_PREFIX_SCAN);
    int gridSize = ceil((float) n / blockSize);
    scan2<add><<<gridSize, blockSize>>>(in, out, sums, n);
    scan2<add><<<1, gridSize>>>(sums, sums2, in, gridSize);
    applySum<add><<<gridSize, blockSize>>>(out, out, sums2, n);
}


//note this works for arrays with up to blocksize^2 elems

void plusScanI(int *in, int *out, int *sums, int *sums2, int n) {
    launchPlusScan(in, out, sums, sums2, n);
}


/*
//note this works for arrays with up to blocksize^2 elems
void plusScanInt(int *in, int *out, int *sums, int *sums2, int n) 
{
   int blockSize = min(n, SHARED_PREFIX_SCAN); 
   int gridSize = ceil ((float) n/blockSize);
   scan2<add><<<gridSize,blockSize>>>(in, out, sums, n);
   scan2<add><<<1, gridSize>>>(sums, sums2, in, gridSize);
   applySum<add><<<gridSize,blockSize>>>(out, out, sums2, n);
}
 */
/*
template <Op op>
void prefixScan(double *in, double *out,
 double *scratch, double *scratch1, int n )
{
   int blockSize;
   int gridSize;
   int gridSizeOld;
   blockSize=min(SHARED_PREFIX_SCAN, n);
   gridSize=(int)ceil((float) n/blockSize);  
   cudaMemset(scratch,0, gridSize*sizeof(double));
   cudaMemset(scratch1,0, gridSize*sizeof(double));
   int i=0;
   plusScan(in, out, scratch, scratch1,n);  
 
}
 */
