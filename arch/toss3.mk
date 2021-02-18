# $Id: makefile.h 1149 2009-12-18 17:35:51Z draeger1 $ 
# Make configuration for toss3 (Livermore)

ARCHDESC = Linux x86-64 AMD Opteron, gcc compilers

#set this variable if an include file is used for extra dependencies
ARCHINCLUDE = 1

CPP = mpicxx
CC = mpicc -std=gnu99
MODULES = "fftw/3.3.8"
WHICH_FFT = USE_FFTW3

# FFT library needs to be included here
CFLAGS_BASE = -D$(WHICH_FFT) -DWITH_MPI -DWITH_PIO -D_GNU_SOURCE \
	      -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas -qopenmp -I /opt/openmpi/3.0/gnu/include/ \
              -I /usr/tce/packages/fftw/fftw-3.3.7-mvapich2-2.2-intel-18.0.1/include/ 
CPPFLAGS_BASE = $(CFLAGS_BASE)
CPPFLAGS += -std=c++0x 

LDFLAGS_BASE = -lm -lc -lstdc++  -Wl,--auto_rpath  -L /opt/openmpi/3.0/gnu/lib/ -l mpi /usr/tce/packages/fftw/fftw-3.3.7-mvapich2-2.2-intel-18.0.1/lib/libfftw3.a 

HAVE_GSL = 0
ifeq ($(HAVE_GSL),1)
CFLAGS_BASE  += -DHAVE_GSL
LDFLAGS_BASE += -lgsl -lgslcblas
endif

HAVE_QHULL = 0
ifeq ($(HAVE_QHULL),1)
  LIBQHULL        = ../base/qhull/libqhull_toss3.a
  QHULLINCLUDEDIR = ../base/qhull/src
  CFLAGS_BASE += -DUSE_QHULL -I$(QHULLINCLUDEDIR)
endif

ifeq ($(HOSTNAME),surface)  
  USE_GPU=1
endif
ifeq ($(HOSTNAME),pascal)  
  USE_GPU=1
  CUDAHOME=/usr/tce/packages/cuda/cuda-10.1.168/
endif

ifeq ($(USE_GPU), 1)
  CFLAGS_BASE += -I $(CUDAHOME)/include/ 
  CFLAGS_BASE += -DUSE_CUDA=1  -DUSE_GPU=1
  LDFLAGS_BASE += -L $(CUDAHOME)/lib64 -lcudart -lcuda -lcurand -lnvrtc -lnvToolsExt -Wl,-rpath,$(CUDAHOME)/lib64

  NVCCFLAGS_BASE = -DUSE_GPU -I ../cub/ --gpu-architecture=sm_60 -I /opt/openmpi/3.0/gnu/include/

  NVCC = $(CUDAHOME)/bin/nvcc -std=c++11

  NVCCFLAGS_OPT  = $(NVCCFLAGS_BASE)  -O3
  NVCCFLAGS_DEBUG  = $(NVCCFLAGS_BASE)  -g -G
  NVCCFLAGS_PROF  = $(NVCCFLAGS_BASE)   -g -pg

endif


CFLAGS_OPT = $(CFLAGS_BASE) -O3 -g
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0 -DLINUX_DEBUG
CFLAGS_PROF = $(CFLAGS_BASE) -DNDEBUG -O3 -g -DPROFILE -pg -fno-inline
CFLAGS_OMP = $(CFLAGS_OPT)

LDFLAGS_OPT = $(LDFLAGS_BASE)
LDFLAGS_DEBUG = $(LDFLAGS_BASE)
LDFLAGS_PROF = $(LDFLAGS_BASE)
LDFLAGS_OMP = $(LDFLAGS_OPT)


ifeq ($(BUILD_MODE),OMP)
  CC += -fopenmp -DWITH_OMP
endif
