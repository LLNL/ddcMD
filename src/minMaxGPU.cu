//#include "minMax.h"
#include <cuda_runtime.h>
#include "cudaUtils.h" 
//#include "thrust/extrema.h"
//#include <thrust/device_vector.h>
#include "system.h"
#include "state.h"
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include "minMaxGPU.h"
#include <utility>
#include <cub/cub.cuh> 
#include <cub/util_allocator.cuh>
#include "gpu_allocator.hpp"
#define CHECK(x) 
#define SHARED_RED 64
typedef double(*Op)(double, double);

//following function from
//http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf
//this fuses all reduction into 1 kernel

__global__ void reduce1XYZ(double *g_idata_xMax, double *g_idata_yMax, double *g_idata_zMax,
                           double *g_idata_xMin, double *g_idata_yMin, double *g_idata_zMin,
                           double *g_odata_xMax, double *g_odata_yMax, double *g_odata_zMax,
                           double *g_odata_xMin, double *g_odata_yMin, double *g_odata_zMin, double * mmgpu, int n)
{
    __shared__ double sdataMax_x[SHARED_RED];
    __shared__ double sdataMax_y[SHARED_RED];
    __shared__ double sdataMax_z[SHARED_RED];
    __shared__ double sdataMin_x[SHARED_RED];
    __shared__ double sdataMin_y[SHARED_RED];
    __shared__ double sdataMin_z[SHARED_RED];

    // each thread loads one element from global to shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
    {
        sdataMax_x[tid] = g_idata_xMax[i];
        sdataMax_y[tid] = g_idata_yMax[i];
        sdataMax_z[tid] = g_idata_zMax[i];

        sdataMin_x[tid] = g_idata_xMin[tid];
        sdataMin_y[tid] = g_idata_yMin[tid];
        sdataMin_z[tid] = g_idata_zMin[tid];
    }
    else
    {
        sdataMax_x[tid] = DBL_MIN;
        sdataMax_y[tid] = DBL_MIN;
        sdataMax_z[tid] = DBL_MIN;

        sdataMin_x[tid] = DBL_MAX;
        sdataMin_y[tid] = DBL_MAX;
        sdataMin_z[tid] = DBL_MAX;
    }
    __syncthreads();

    // do reduction in shared mem
    for (unsigned int s = 1; s < blockDim.x; s *= 2)
    {
        if (tid % (2 * s) == 0)
        {
            sdataMin_x[tid] = fmin(sdataMin_x[tid], sdataMin_x[tid + s]);
            sdataMin_y[tid] = fmin(sdataMin_y[tid], sdataMin_y[tid + s]);
            sdataMin_z[tid] = fmin(sdataMin_z[tid], sdataMin_z[tid + s]);
            sdataMax_x[tid] = fmax(sdataMax_x[tid], sdataMax_x[tid + s]);
            sdataMax_y[tid] = fmax(sdataMax_y[tid], sdataMax_y[tid + s]);
            sdataMax_z[tid] = fmax(sdataMax_z[tid], sdataMax_z[tid + s]);
        }
        __syncthreads();
    }
    // write result for this block to global mem
    if (tid == 0)
    {
        g_odata_xMin[blockIdx.x] = sdataMin_x[0];
        g_odata_yMin[blockIdx.x] = sdataMin_y[0];
        g_odata_zMin[blockIdx.x] = sdataMin_z[0];
        g_odata_xMax[blockIdx.x] = sdataMax_x[0];
        g_odata_yMax[blockIdx.x] = sdataMax_y[0];
        g_odata_zMax[blockIdx.x] = sdataMax_z[0];

        mmgpu[0] = sdataMin_x[0];
        mmgpu[1] = sdataMax_x[0];
        mmgpu[2] = sdataMin_y[0];
        mmgpu[3] = sdataMax_y[0];
        mmgpu[4] = sdataMin_z[0];
        mmgpu[5] = sdataMax_z[0];
    }
}

template <Op op> __global__ void reduce1(double *g_idata, double *g_odata, int n)
{
    __shared__ double sdata[SHARED_RED];
    // each thread loads one element from global to shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x
        * blockDim.x + threadIdx.x;
    if (i >= n)
    {
        return;
    }
    sdata[tid] = g_idata[i];
    __syncthreads();
    // do reduction in shared mem
    for (unsigned int s = 1; s < blockDim.x; s *= 2)
    {
        if (tid % (2 * s) == 0)
        {
            sdata[tid] = op(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

void execReduceXYZ(double *g_idata_x, double *g_idata_y, double *g_idata_z,
                   double *g_odata_xMax, double *g_odata_yMax, double *g_odata_zMax,
                   double *g_odata_xMin, double *g_odata_yMin, double *g_odata_zMin, double *mmgpu, int n)
{
    int blockSize = SHARED_RED;
    int gridSize = (int) ceil((float) n / blockSize);
    //int gridSizeOld;
    reduce1XYZ << <gridSize, blockSize>>>(g_idata_x, g_idata_y, g_idata_z,
        g_idata_x, g_idata_y, g_idata_z,
        g_odata_xMax, g_odata_yMax, g_odata_zMax,
        g_odata_xMin, g_odata_yMin, g_odata_zMin, mmgpu, n);

    while (gridSize > 1)
    {
        //gridSizeOld=gridSize;
        gridSize = (int) ceil((float) gridSize / blockSize);
        //printf("gridsize %i\n", gridSize);
        reduce1XYZ << <gridSize, blockSize>>>(g_odata_xMax, g_odata_yMax, g_odata_zMax,
            g_odata_xMin, g_odata_yMin, g_odata_zMin,
            g_odata_xMax, g_odata_yMax, g_odata_zMax,
            g_odata_xMin, g_odata_yMin, g_odata_zMin, mmgpu, n);
    }
}

template <Op op>
void execReduce(double *array, STATE* hostState, double *mmgpu, int offset,
                double *scratch, double *scratch1)
{
    int blockSize;
    int gridSize;
    int gridSizeOld;
    blockSize = min(SHARED_RED, hostState->nlocal);
    gridSize = (int) ceil((float) hostState->nlocal / blockSize);
    CHECK(printf("blocksize %i gridsize %i\n", blockSize, gridSize);)
    cudaMemset(scratch, 0, gridSize * sizeof (double));
    cudaMemset(scratch1, 0, gridSize * sizeof (double));
    int i = 0;

    //each iteration stores the min/max element from each block "b" into 
    //output scratch space array at position array[b]. 
    //the output scratch space array then becomes the input array
    //in the next iteration. this is implemented via a pointer swap
    //the last iteration will thus be that in which the gridsize=1 (only 1 block)
    reduce1<op> << <gridSize, blockSize>>>(array, scratch, hostState->nlocal);
    while (gridSize > 1)
    {
        blockSize = min(SHARED_RED, gridSize);
        gridSizeOld = gridSize;
        gridSize = (int) ceil((float) gridSize / blockSize);
        CHECK(printf("blocksize %i gridsize %i\n", blockSize, gridSize);)

        reduce1<op> << <gridSize, blockSize>>>(scratch, scratch1, gridSizeOld);
        std::swap(scratch, scratch1);
        i++;
    }
    double *mmgpuOffset = mmgpu + offset;
    gpu_memcpy_device2device(mmgpuOffset, scratch, 1);
    if (i % 2)
    {
        std::swap(scratch, scratch1);
    }
}

/*
double *getMinMaxGPU(COLLECTION *collection, STATE *gpuState_h, STATE* hostState, double *mmgpu)
{
   execReduceXYZ(gpuState_h->rx, gpuState_h->ry, gpuState_h->rz,
                 collection->scratch, collection->scratch1, collection->scratch2,
//                 collection->scratch3, collection->scratch4, collection->scratch5, mmgpu, hostState->nion);
   
}*/

double *getMinMaxGPU(COLLECTION *collection, STATE *gpuState_h, STATE* hostState, double *mmgpu) 
{
    //GPUNLIST *gnlist = collection->gnlist; 
    //cub::CachingDeviceAllocator  g_allocator(true);
    size_t temp_storage_bytes = 0;
    void *d_temp_storage = NULL;
    //size_t avail_bytes = gnlist->allocFactor*hostState->nion*sizeof(double);
    cub::DeviceReduce::Min(d_temp_storage, temp_storage_bytes, gpuState_h->rx, mmgpu, hostState->nion);
    //g_allocator.DeviceAllocate(&d_temp_storage, temp_storage_bytes);	
    cudaMalloc(&d_temp_storage, temp_storage_bytes);

    cub::DeviceReduce::Min(d_temp_storage, temp_storage_bytes, gpuState_h->rx, mmgpu, hostState->nion);
    cub::DeviceReduce::Max(d_temp_storage, temp_storage_bytes, gpuState_h->rx, mmgpu + 1, hostState->nion);
    cub::DeviceReduce::Min(d_temp_storage, temp_storage_bytes, gpuState_h->ry, mmgpu + 2, hostState->nion);
    cub::DeviceReduce::Max(d_temp_storage, temp_storage_bytes, gpuState_h->ry, mmgpu + 3, hostState->nion);
    cub::DeviceReduce::Min(d_temp_storage, temp_storage_bytes, gpuState_h->rz, mmgpu + 4, hostState->nion);
    cub::DeviceReduce::Max(d_temp_storage, temp_storage_bytes, gpuState_h->rz, mmgpu + 5, hostState->nion);

    //if (d_temp_storage) g_allocator.DeviceFree(d_temp_storage);
    if (d_temp_storage) cudaFree(d_temp_storage);
    CUDA_SAFE_CALL(cudaPeekAtLastError());

    return mmgpu;
}
//this ran 6 sepearte reduce kernels to find mins and maxs
//in x y and z. it results in to much kernel launch/driver overhead

double* getMinMaxGPUOld(COLLECTION *collection, STATE * gpuState_h, STATE* hostState, double *mmgpu) 
{
     GPUNLIST *gnlist = collection->gnlist; 
    double *scratch = gnlist->scratch;
    double *scratch1 = gnlist->scratch1;
    //int blockSize, gridSize;
    //blockSize=min(SHARED_RED, hostState->nion);
    //gridSize=(int)ceil((float)hostState->nion/blockSize);

    //allocate scratch space for reduction kernels

    //execute reduction kernel routines
    execReduce<fmin>(gpuState_h->rx, hostState, mmgpu, 0, scratch, scratch1);
    execReduce<fmax>(gpuState_h->rx, hostState, mmgpu, 1, scratch, scratch1);
    execReduce<fmin>(gpuState_h->ry, hostState, mmgpu, 2, scratch, scratch1);
    execReduce<fmax>(gpuState_h->ry, hostState, mmgpu, 3, scratch, scratch1);
    execReduce<fmin>(gpuState_h->rz, hostState, mmgpu, 4, scratch, scratch1);
    execReduce<fmax>(gpuState_h->rz, hostState, mmgpu, 5, scratch, scratch1);

    //int h = hostState->nlocal;
    //double lx, ly, lz, hx,hy,hz;
    CHECK(
          for (int i = 0; i < h; i++)
    {
          vlx = std::min(lx, hostState->rx[i]);
          hx = std::max(hx, hostState->rx[i]);
          ly = std::min(ly, hostState->ry[i]);
          hy = std::max(hy, hostState->ry[i]);
          lz = std::min(lz, hostState->rz[i]);
          hz = std::max(hz, hostState->rz[i]);

    }
          printf("lx  %f hx %f ly %f hy %f lz %f hz %f\n", lx, hx, ly, hy, lz, hz);
          )
        //cudaFree(scratch);
        //cudaFree(scratch1);
    return mmgpu;
}
/*
long timeMinMaxThrust(int len)
{
    struct timeval start, end;
    thrust::device_vector<double> vec(len);
    for (int ijk =0; ijk <len; ijk++)
    {   
       vec[ijk] = ijk;    
    }   
    //CuAssertTrue(tc,1==1);
    

    thrust::device_vector<double> d_vec = vec;
    gettimeofday(&start, NULL);
    thrust::device_vector<double>::iterator iter = thrust::max_element(d_vec.begin(), d_vec.end());
    gettimeofday(&end, NULL);
    unsigned int position = iter - d_vec.begin();
    double max_val = *iter;

    long spent = (end.tv_sec * 1000000 + end.tv_usec) -(start.tv_sec * 1000000 + start.tv_usec);
    printf("thrust time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
                  - (start.tv_sec * 1000000 + start.tv_usec)));
    std::cout << "The maximum value is " << max_val << " at position " << position << std::endl;
    
    // empty the vector
    d_vec.clear();

    // deallocate any capacity which may currently be associated with vec
    d_vec.shrink_to_fit();

    return spent;
}*/
