#ifndef GPU_ALLOCATOR_HPP
#define GPU_ALLOCATOR_HPP

#include <cuda_runtime_api.h>
#include <cuda.h>

#include <execinfo.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>

static void full_write(int fd, const char *buf, size_t len)
{
        while (len > 0) {
                ssize_t ret = write(fd, buf, len);

                if ((ret == -1) && (errno != EINTR))
                        break;

                buf += (size_t) ret;
                len -= (size_t) ret;
        }
}

inline void print_backtrace(void)
{
        static const char start[] = "BACKTRACE ------------\n";
        static const char end[] = "----------------------\n";

        void *bt[1024];
        int bt_size;
        char **bt_syms;
        int i;

        bt_size = backtrace(bt, 1024);
        bt_syms = backtrace_symbols(bt, bt_size);
        full_write(STDERR_FILENO, start, strlen(start));
        for (i = 1; i < bt_size; i++) {
                size_t len = strlen(bt_syms[i]);
                full_write(STDERR_FILENO, bt_syms[i], len);
                full_write(STDERR_FILENO, "\n", 1);
        }
        full_write(STDERR_FILENO, end, strlen(end));
    free(bt_syms);
}

template<typename T> void gpu_allocator(T &ptrref, long long int n) {
  //print_backtrace();
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
