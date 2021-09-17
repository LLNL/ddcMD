## INSTALLATION GUIDE 

## 1. Obtain the code
Download the code from GitHub (to make git submodule work user has to do the git clone to download the code)
```
git clone git@github.com:LLNL/ddcMD.git
```

## 2. Build code
After downloading, go to the ddcMD home directory and use git submodule command to obtain three dependent libraries: 
NVIDIA [CUB](https://github.com/NVlabs/cub) library and LLNL [simutil](https://github.com/LLNL/simutil) and [recbis](https://github.com/LLNL/recbis) libraries.
```
cd ddcMD
git submodule update --init --recursive
```

### 2.1. Use CMake 
Compile code without GPU/CUDA

```
mkdir build && cd build
cmake -DUSE_GPU=OFF ../
make -j16
```

Compile code with GPU/CUDA

```
mkdir build && cd build
cmake -DUSE_GPU=ON ../
make -j16
```

Use IBM XL compiler
```
CC=xlc CXX=xlc++ cmake ../
```

Use Intel compiler

```
CC=icc CXX=icpc cmake ../
```

### 2.2. Use makefile

The makefile is in ddcMD/src directory. If the make command is successful, an executable, ddcmd-{arch}, will be generated in ddcMD/bin
```
cd src
make
```

#### 2.2.1. Support architectures
Currently, the code only supports a few architectures as shown in ddcMD/arch
```
apple.mk  armbuntu.mk  macosx.mk  sierra.mk  summit.mk  toss3.mk
```

The code determines the architecture by the hostname of the machine defined in ddcMD/arch/Makefile.arch. For example, in the ddcMD/arch/Makefile.arch the script to build code on the Summit at ORNL is defined as
```
ifeq ($(HOSTNAME_D), summit)
  ARCHGUESS = summit
endif
```
The script will invoke summit.mk during the building process.

If user wants to build the code on a machine with the hostname that cannot find in ddcMD/arch/Makefile.arch, user can add the new hostname to the file and use the corresponding \*.mk architecture file. If user cannot find the architecture file in ddcMD/arch, user can use any of the \*.mk files as template and modify compiler flags and libraries.


