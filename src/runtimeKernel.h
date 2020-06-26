#ifndef RUNTIME_KERNEL_H
#define RUNTIME_KERNEL_H
#include <string>
#include <nvrtc.h>
#include <cuda.h>
#include <vector>
#include <map>

typedef std::pair<std::vector<std::string>, std::map<std::string, int> > runtimeKernelSegments;

class RuntimeKernel
{
public:
    nvrtcProgram prog;
    char * ptx;
    nvrtcResult compileResult;
    std::string name; //name given to kernel
    std::string kernelBody; //body of generated kernel
    CUmodule module; //module into which ptx of new kernel will be loaded
};

class RuntimeKernelManager
{
    CUdevice cuDevice;
    CUcontext context;

    std::vector<RuntimeKernel> kernels;
public:
    void compile(std::string fileName, std::string kernelName);
    void compileFromString(std::string kernelBody, std::string kernelName);
    void loadToRuntime(std::string kernelName);
    void getKernel(CUfunction *func, std::string kernelName);
    void initDeviceContext();
    runtimeKernelSegments parseKernelSegmentFile(std::string templateFileName);
    std::string generateKernel(runtimeKernelSegments &segments, void * parms);
};

#endif
