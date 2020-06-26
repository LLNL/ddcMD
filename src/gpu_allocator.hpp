#ifndef GPU_ALLOCATOR_HPP
#define GPU_ALLOCATOR_HPP

#include <cuda_runtime_api.h>
#include <cuda.h>

template<typename T> void gpu_allocator(T &ptrref, long long int n) {
  CUDA_SAFE_CALL( cudaMalloc((void **) &ptrref, sizeof(*ptrref) * n);)
}

template<typename T1, typename T2> void gpu_memcpy_host2device(T1 &ptrref,  T2 &srcref, long long int n) {
  CUDA_SAFE_CALL(cudaMemcpy(ptrref, srcref, sizeof(*ptrref)*n, cudaMemcpyHostToDevice);)
}

template<typename T1, typename T2> void gpu_memcpy_device2host(T1 &ptrref,  T2 &srcref, long long int n) {
  CUDA_SAFE_CALL(cudaMemcpy(ptrref, srcref, sizeof(*ptrref)*n, cudaMemcpyDeviceToHost);)
}

template<typename T1, typename T2> void gpu_memcpy_device2device(T1 &ptrref,  T2 &srcref, long long int n) {
  CUDA_SAFE_CALL(cudaMemcpy(ptrref, srcref, sizeof(*ptrref)*n, cudaMemcpyDeviceToDevice);)
}

#endif //GPU_ALLOCATOR_HPP
