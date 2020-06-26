//generates fast constraint kernels
#include "runtimeKernel.h"
#include <iostream>
#include <fstream>
#include <string>
#include <nvrtc.h>
#include <cuda.h>
#include <streambuf>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include "bioCharmm.h"
#define NVRTC_SAFE_CALL(x)                                        \
  do {                                                            \
    nvrtcResult result = x;                                       \
    if (result != NVRTC_SUCCESS) {                                \
      std::cerr << "\nerror: " #x " failed with error "           \
                << nvrtcGetErrorString(result) << '\n';           \
      exit(1);                                                    \
    }                                                             \
  } while(0)
#define CUDA_SAFE_CALL2(x)                                         \
  do {                                                            \
    CUresult result = x;                                          \
    if (result != CUDA_SUCCESS) {                                 \
      const char *msg;                                            \
      cuGetErrorName(result, &msg);                               \
      std::cerr << "\nerror: " #x " failed with error "           \
                << msg << '\n';                                   \
      exit(1);                                                    \
    }                                                             \
  } while(0)

void RuntimeKernelManager::initDeviceContext()
{

    CUDA_SAFE_CALL2(cuInit(0));
    CUDA_SAFE_CALL2(cuDeviceGet(&cuDevice, 0));
    CUDA_SAFE_CALL2(cuCtxCreate(&context, 0, cuDevice));

}

void RuntimeKernelManager::compile(std::string fileName, std::string kernelName)
{
    printf("compiling \n");
    std::ifstream t(fileName.c_str());
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

    nvrtcProgram p;
    // Create an instance of nvrtcProgram with the SAXPY code string.
    NVRTC_SAFE_CALL(
                    nvrtcCreateProgram(&p, // prog
                                       str.c_str(), // buffer
                                       fileName.c_str(), // name
                                       0, // numHeaders
                                       NULL, // headers
                                       NULL)); // includeNames
    const char *opts[] = {"--gpu-architecture=compute_60",
        "-default-device"};
    nvrtcResult compileResult = nvrtcCompileProgram(p, // prog
                                                    2, // numOptions
                                                    opts); // options
    RuntimeKernel k;
    k.prog = p; //std::move(p);
    k.compileResult = compileResult;
    k.kernelBody = str;
    k.name = kernelName;
    kernels.push_back(k); //std::move(k));

    // Obtain compilation log from the program.
    size_t logSize;
    NVRTC_SAFE_CALL(nvrtcGetProgramLogSize(p, &logSize));
    char *log = new char[logSize];
    NVRTC_SAFE_CALL(nvrtcGetProgramLog(p, log));
    std::cout << log << '\n';
    delete[] log;
    if (compileResult != NVRTC_SUCCESS)
    {
        printf("bad compile for %s\n", fileName.c_str());
        exit(1);
    }
}

void RuntimeKernelManager::loadToRuntime(std::string kernelName)
{
    int kernelId = 0;
    for (int i = 0; i < kernels.size(); i++)
    {
        if (kernels[i].name == kernelName)
        {
            kernelId = i;
            printf("found kernel \n");
        }
    }

    RuntimeKernel *k = &kernels[kernelId];

    // Obtain PTX from the program.
    size_t ptxSize;
    NVRTC_SAFE_CALL(nvrtcGetPTXSize(k->prog, &ptxSize));
    k->ptx = new char[ptxSize];
    NVRTC_SAFE_CALL(nvrtcGetPTX(k->prog, k->ptx));
    //we don't need this anymore
    NVRTC_SAFE_CALL(nvrtcDestroyProgram(&k->prog));
    CUDA_SAFE_CALL2(cuModuleLoadDataEx(&k->module, k->ptx, 0, 0, 0));
}

void RuntimeKernelManager::getKernel(CUfunction * func, std::string kernelName)
{
    int kernelId = 0;
    for (int i = 0; i < kernels.size(); i++)
    {
        if (kernels[i].name == kernelName)
        {
            kernelId = i;
            printf("found kernel \n");
        }
    }


    CUDA_SAFE_CALL2(cuModuleGetFunction(func, kernels[kernelId].module, kernelName.c_str()));
}

runtimeKernelSegments RuntimeKernelManager::parseKernelSegmentFile(std::string fileName)
{
    printf("parsing \n");
    std::ifstream testFile(fileName);
    std::string currentSection = "";
    std::string codeSegment = "";
    int sectionCounter = 0;

    std::vector<std::string> sections; //strings that compose kernel body, modified by generateKernel()
    std::map<std::string, int> sectionIDMap; //map of section name to section id in "sections"
    std::string line;
    std::string delim = " ";
    //loop through each line of file
    //and find all listed code "segments", where each segment in file is 
    //the continguous chunk of lines between pair of "@startJit [SECTIONNAME]" and "@endJit [SECTIONNAME]" keywords
    //then load each section into a vector
    //also keep a hashmap that maps each segment's SECTIONNAME to it's index in the vector
    //just so it's more straightforward to access code segments
    while (std::getline(testFile, line))
    {
        std::istringstream split(line);
        std::vector<std::string> tokens;
        std::string item;
        while (std::getline(split, item, ' '))
        {
            tokens.push_back(item);
        }
        if (tokens.size() > 0)
        {
            if (tokens[0] == "@startJit")
            {
                //get name of current code segment, insert into code segment id->name map
                currentSection = tokens[1];
                printf("parsing %s\n", currentSection.c_str());
                sectionIDMap.insert(std::pair<std::string, int>(currentSection, sectionCounter));
                sectionCounter++;
            }
            else if (tokens[0] == "@endJit")
            {
                //insert code segment to our code segment map, reset codesegment buffer
                sections.push_back(codeSegment + "\n");
                codeSegment = "";
                printf("finished parsing %s\n", currentSection.c_str());
            }
            else if (currentSection != "")
            {
                //add current line to code segment buffer
                codeSegment += line + "\n";
            }
        }
    }

    int loopID = sectionIDMap["loop"];
    //std::pair<std::vector<std::string>, std::map<std::string, int>> runtimeKernelSegments;
    return std::pair<std::vector<std::string>, std::map < std::string, int>>(sections, sectionIDMap);
    //printf("section loop %s\n", sections[loopID].c_str());
}

std::string RuntimeKernelManager::generateKernel(runtimeKernelSegments &segments, void *parms)
{
    /*
    switch(resTypeID)
   {
        for (;it<maxit;it++) //start the iterative loop
        {
            errMax =0.0;
            for (int ab=0; ab<numPair; ab++)
            {
                int pairID = pairIDs[ab];
                int a = pairsIindexs[ab];
                int b = pairsJindexs[ab];
                double dist2 = distances2[ab];
                double dist2inv = distances2inv[ab];
                THREE_VECTOR vab; 
                VOP2(vab,=,v[a],-,v[b]);
                double rma=rMass[a];
                double rmb=rMass[b];                           
                double rma_rmb_inv = rMassPairInv[ab];
                //call Front or back function
                double rvab = FrontFuncOpt(dt, dt2inv,dist2, rab[ab], vab)*dist2inv;
                //if (wid==0) printf("it %i ab %idist2 %f rvab %f dt %f nc %i\n", it, ab,dist2, rvab, dt, nConstraint);
                double gab=-rvab*rma_rmb_inv;   //units: mass/time)
                double err = fabs(rvab*dt);
                errMax = fmax(err, errMax);
                VSVOP(v[a],+=,(rma*gab),*,rab[ab]); 
                VSVOP(v[b],-=,(rmb*gab),*,rab[ab]); 
                gamma[ab] += gab;  
            } //numPair
            if(tid==0 && errMax<tol) {break;}
        } //maxit
   //      break;
  // }
   }

     */
    CHARMMPOT_PARMS *cp_parms = static_cast<CHARMMPOT_PARMS*> (parms);
    CHARMM_PARMS *charmmParms = cp_parms->charmmParms;
    std::string tab = "    ";
    std::string curtab = "";
    std::ostringstream curr;
    curr << curtab + "if(tid==0){\n";
    std::cout << std::fixed;
    int precision = 17;
    curtab = curtab + tab;
    curr << curtab + "switch(resTypeID){\n";
    curtab = curtab + tab;

    //int constraintListSize = 0; //also a counter for making a prefix sum array
    for (int i = 0; i < charmmParms->resiConnSize; i++)
    {
        if (i != 2)continue;
        curr << curtab + "case " + std::to_string(i) + ":\n";
        RESI_CONN * resiConn = charmmParms->resiConnList[i];
        CONSTRAINT ** consList = resiConn->consList;
        //iterate over constraint groups in a residue
        if (resiConn->consListSize > 1)continue;
        for (int j = 0; j < resiConn->consListSize; j++)
        {
            CONSTRAINT *constraint = consList[j];
            //iterate over pairs in a constraint group
            if (constraint->numPair > 0)
            {
                printf("i %i j %i\n", i, j);
                curtab = curtab + tab;
                //curr<< curtab+"int it=0;\n";
                curr << curtab + "for (int it=0;it<maxit;it++){\n";
                curtab = curtab + tab;
                curr << curtab + "errMax =0.0;\n";
                curr << curtab + "int ab=0;\n";
                for (int k = 0; k < constraint->numPair; k++)
                {
                    curr << curtab + "{\n";
                    curtab = curtab + tab;
                    CONS_PAIR *cpair = constraint->conspairList[k];
                    double dist = cpair->distance;
                    int a = cpair->atomIindex;
                    int b = cpair->atomJindex;
                    double mass_a = resiConn->atomList[a]->atmTypePtr->mass;
                    double mass_b = resiConn->atomList[b]->atmTypePtr->mass;
                    curr << curtab + "const int a = " << std::to_string(a) + ";\n";
                    curr << curtab + "const int b = " + std::to_string(b) + ";\n";
                    curr << curtab + "const double dist2 = " << std::setprecision(precision) << dist * dist << ";\n";
                    curr << curtab + "const double dist2inv = " << std::setprecision(precision) << 1.0 / (dist * dist) << ";\n";
                    curr << curtab + "const double rma = " << std::setprecision(precision) << mass_a << ";\n";
                    curr << curtab + "const double rmb = " << std::setprecision(precision) << mass_b << ";\n";
                    curr << curtab + "const double rma_rmb_inv = " << std::setprecision(precision) << 1.0 / (mass_a + mass_b) << ";\n";
                    curr << curtab + "THREE_VECTOR vab;\n";
                    curr << curtab + "VOP2(vab,=,v[a],-,v[b]);\n";
                    curr << curtab + "double rvab = FrontFuncOpt(dt, dt2inv,dist2, rab[ab], vab)*dist2inv;\n";
                    curr << curtab + "double gab=-rvab*rma_rmb_inv;\n";
                    curr << curtab + "double err = fabs(rvab*dt);\n";
                    curr << curtab + "errMax = fmax(err, errMax);\n";
                    curr << curtab + "VSVOP(v[a],+=,(rma*gab),*,rab[ab]);\n";
                    curr << curtab + "VSVOP(v[b],-=,(rmb*gab),*,rab[ab]);\n";
                    curr << curtab + "gamma[ab] += gab;\n";
                    curtab = curtab.substr(0, curtab.length() - 4);
                    curr << curtab + "}\n";
                    curr << curtab + "ab++;\n";
                }
                curr << curtab + "if(tid==0 && errMax<tol) {break;}\n";
                curtab = curtab.substr(0, curtab.length() - 2 * 4);
                curr << curtab + "};\n";
            }
        }
        curr << curtab + "break; //maxit \n";
    }
    //defaulta
    curr << curtab + "default:\n";
    std::string defaultcase = "\
             for (;it<maxit;it++) //start the iterative loop\n\
             {\n\
                errMax =0.0;\n\
                for (int ab=0; ab<numPair; ab++)\n\
                {\n\
                   int pairID = pairIDs[ab];\n\
                   int a = pairsIindexs[ab];\n\
                   int b = pairsJindexs[ab];\n\
                   double dist2 = distances2[ab];\n\
                   double dist2inv = distances2inv[ab];\n\
                   THREE_VECTOR vab; \n\
                   VOP2(vab,=,v[a],-,v[b]);\n\
                   double rma=rMass[a];\n\
                   double rmb=rMass[b]; \n\
                   double rma_rmb_inv = rMassPairInv[ab];\n\
                   //call Front or back function\n\
                   double rvab = FrontFuncOpt(dt, dt2inv,dist2, rab[ab], vab)*dist2inv;\n\
                   double gab=-rvab*rma_rmb_inv;   //units: mass/time)\n\
                   double err = fabs(rvab*dt);\n\
                   errMax = fmax(err, errMax);\n\
                   VSVOP(v[a],+=,(rma*gab),*,rab[ab]); \n\
                   VSVOP(v[b],-=,(rmb*gab),*,rab[ab]); \n\
                   gamma[ab] += gab;  \n\
                } //numPair\n\
                if(tid==0 && errMax<tol) {break;}\n\
             } //maxit\n";
    curr << defaultcase;
    curtab = curtab.substr(0, curtab.length() - 4);
    curr << curtab + "}\n";
    curtab = curtab.substr(0, curtab.length() - 4);
    curr << curtab + "}\n";
    printf("done iter\n");
    //printf("curr \n%s\n", curr.c_str());
    auto codeSegments = segments.first;
    auto sectionIDMap = segments.second;
    int loopID = sectionIDMap["loop"];
    codeSegments[loopID] = curr.str();
    std::string fileString;
    for (auto segment : codeSegments)
    {
        std::cout << segment << "\n";
        fileString += segment;
    }
    std::ofstream out("outputKernel.cu");
    out << fileString;
    out.close();
    return fileString;
}

int main1()
{
    RuntimeKernelManager kernelGen;
    kernelGen.compile("kernelTemplate.cu", "");
    //init the runtime
    kernelGen.initDeviceContext();
    // Load the generated PTX and get a handle to the SAXPY kernel.
    kernelGen.loadToRuntime("kernelTemplate.cu");

    //only this should be called in the simulation loop
    //CUfunction kernel;
    //kernelGen.getKernel(kernel, "saxpy");

    /*
    // Generate input for execution, and create output buffers.
    size_t n = NUM_THREADS * NUM_BLOCKS;
    size_t bufferSize = n * sizeof(float);
    float a = 5.1f;
    float *hX = new float[n], *hY = new float[n], *hOut = new float[n];
    for (size_t i = 0; i < n; ++i) {
      hX[i] = static_cast<float>(i);
      hY[i] = static_cast<float>(i * 2);
    }
    CUdeviceptr dX, dY, dOut;
    CUDA_SAFE_CALL(cuMemAlloc(&dX, bufferSize));
    CUDA_SAFE_CALL(cuMemAlloc(&dY, bufferSize));
    CUDA_SAFE_CALL(cuMemAlloc(&dOut, bufferSize));
    CUDA_SAFE_CALL(cuMemcpyHtoD(dX, hX, bufferSize));
    CUDA_SAFE_CALL(cuMemcpyHtoD(dY, hY, bufferSize));
    // Execute SAXPY.
    void *args[] = { &a, &dX, &dY, &dOut, &n };
    CUDA_SAFE_CALL(
      cuLaunchKernel(kernel,
                     NUM_BLOCKS, 1, 1,    // grid dim
                     NUM_THREADS, 1, 1,   // block dim
                     0, NULL,             // shared mem and stream
                     args, 0));           // arguments
    CUDA_SAFE_CALL(cuCtxSynchronize());
    // Retrieve and print output.
    CUDA_SAFE_CALL(cuMemcpyDtoH(hOut, dOut, bufferSize));
    for (size_t i = 0; i < n; ++i) {
      std::cout << a << " * " << hX[i] << " + " << hY[i]
                << " = " << hOut[i] << '\n';
    }
    // Release resources.
    CUDA_SAFE_CALL(cuMemFree(dX));
    CUDA_SAFE_CALL(cuMemFree(dY));
    CUDA_SAFE_CALL(cuMemFree(dOut));
    //CUDA_SAFE_CALL(cuModuleUnload(module));
    //CUDA_SAFE_CALL(cuCtxDestroy(context));
    delete[] hX;
    delete[] hY;
    delete[] hOut;
     */
    return 0;


}
