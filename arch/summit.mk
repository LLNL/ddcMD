# $Id: makefile.h 1149 2009-12-18 17:35:51Z draeger1 $ 
# Make configuration for Power clusters (Livermore)

ARCHDESC = Linux x86-64 IBM Power + NVIDIA Volta

#set this variable if an include file is used for extra dependencies
ARCHINCLUDE = 1
#GNU options
#CC =/usr/tcetmp/bin/mpigcc -std=gnu99
#CPP = /usr/tcetmp/bin/mpig++ -std=c++11
#FC = /usr/tcetmp/bin/mpigfortran
#XL options
CC = mpixlc
CPP = mpixlC -std=c++11 -qlanglvl=extended0x
FC = mpixlf

CPPL = $(shell dirname `which xlc++`)/../lib/
MPI_INC = $(shell dirname `which mpixlC`)/../include/
MPI_LIB = $(shell dirname `which mpixlC`)/../lib/
#CPPL = /usr/tce/packages/xl/xl-beta-2018.04.06/lib/
#CPPL = /sw/summit/xl/16.1.1-4/xlC/16.1.1/lib/

# FFT library needs to be included here
#-DUSE_FFTW3 
#/usr/tcetmp/packages/fftw/fftw-3.3.5-xl-15.1.5/include
#
CFLAGS_BASE =  -DWITH_MPI -DWITH_PIO -D_GNU_SOURCE -Wall -Wextra -Wno-unused-parameter \
              -Wno-unknown-pragmas
#GNU Compiler
FFTW_DIR=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-5/fftw-3.3.8-azzdjlzx2j6dpqvzdir2nwvxypohyfq4/lib/
LDFLAGS_BASE =  -lstdc++ -lm -lc $(FFTW_DIR)/libfftw3.a -L $(CPPL) -L$(MPI_LIB)
#XL Compiler
#LDFLAGS_BASE =  -libmc++ -lstdc++ -lm -lc /usr/tcetmp/packages/fftw/fftw-3.3.5-xl-15.1.5/lib/libfftw3.a -L $(CPPL)


HAVE_GSL = 0
ifeq ($(HAVE_GSL),1)
CFLAGS_BASE  += -DHAVE_GSL
LDFLAGS_BASE += -lgsl -lgslcblas
endif

HAVE_QHULL = 0
ifeq ($(HAVE_QHULL),1)
  LIBQHULL        = ../base/qhull/libqhull_power.a
  QHULLINCLUDEDIR = ../base/qhull/src
  CFLAGS_BASE += -DUSE_QHULL -I$(QHULLINCLUDEDIR)
endif

USE_GPU = 1
ifeq ($(USE_GPU),1)
  CUDAHOME=/sw/summit/cuda/10.1.168/
  CFLAGS_BASE += -I $(CUDAHOME)/include/ 
  CFLAGS_BASE +=  -DUSE_GPU=1
  LDFLAGS_BASE += -L $(CUDAHOME)/lib64 -lcudart -lcuda -lcurand -lnvToolsExt -lnvrtc -Wl,-rpath,$(CUDAHOME)/lib64

  NVCCFLAGS_BASE += -DUSE_GPU=1 -I ../cub --gpu-architecture=sm_70 -I/usr/local/tools/mvapich-gnu/include -I$(MPI_INC)
  NVCC = $(CUDAHOME)/bin/nvcc  -std=c++11
#    NVCCFLAGS += --gpu-architecture=sm_60 -I/usr/local/tools/mvapich-gnu/include -I /usr/tce/packages/spectrum-mpi/ibm/spectrum-mpi-2018.02.05/include/
endif


CFLAGS_OPT = $(CFLAGS_BASE)   -O3 -g 
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0 -DLINUX_DEBUG
CFLAGS_PROF = $(CFLAGS_BASE) -DNDEBUG -O3 -g -DPROFILE -pg -fno-inline
CFLAGS_OMP = $(CFLAGS_OPT)

FFLAGS_OPT = -O3 -g
FFLAGS_DEBUG = -g -DLINUX_DEBUG -pg -cpp -w 
FFLAGS_PROF = -O3 -g
FFLAGS_OMP = $(FFLAGS_OPT)

NVCCFLAGS_OPT  = $(NVCCFLAGS_BASE)  -O3 
NVCCFLAGS_DEBUG  = $(NVCCFLAGS_BASE)  -g -G 
NVCCFLAGS_PROF  = $(NVCCFLAGS_BASE)   -g -pg 

LDFLAGS_OPT = $(LDFLAGS_BASE)
LDFLAGS_DEBUG = $(LDFLAGS_BASE)
LDFLAGS_PROF = $(LDFLAGS_BASE)
LDFLAGS_OMP = $(LDFLAGS_OPT)


ifeq ($(BUILD_MODE),OMP)
  CC += -fopenmp -DWITH_OMP
endif
