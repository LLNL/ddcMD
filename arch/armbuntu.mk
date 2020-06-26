# $Id: makefile.h 1149 2009-12-18 17:35:51Z draeger1 $ 
# Make configuration for Peloton clusters (Livermore)

ARCHDESC = Linux x86-64 AMD Opteron, gcc compilers

#set this variable if an include file is used for extra dependencies
ARCHINCLUDE = 1

CC = mpicc -std=gnu99
FC = mpif90

# FFT library needs to be included here
#CFLAGS_BASE = -DUSE_FFTW2 -I/usr/local/tools/fftw-2.1.5/include \
#              -DWITH_MPI -DWITH_PIO -Wall -Wextra -Wno-unused-parameter \
#              -Wno-unknown-pragmas
#LDFLAGS_BASE = -lm -lc /usr/local/tools/fftw-2.1.5/lib/libfftw.a 

CFLAGS_BASE = -DUSE_FFTW2 \
              -DWITH_MPI -DWITH_PIO -D_GNU_SOURCE -Wall -Wextra -Wno-unused-parameter \
              -Wno-unknown-pragmas
LDFLAGS_BASE = -lm -lfftw

HAVE_GSL = 0
ifeq ($(HAVE_GSL),1)
CFLAGS_BASE  += -DHAVE_GSL
LDFLAGS_BASE += -lgsl -lgslcblas
endif

HAVE_QHULL = 0
ifeq ($(HAVE_QHULL),1)
  LIBQHULL        = ../base/qhull/libqhull_armbuntu.a
  QHULLINCLUDEDIR = ../base/qhull/src
  CFLAGS_BASE += -DUSE_QHULL -I$(QHULLINCLUDEDIR)
endif

CFLAGS_OPT = $(CFLAGS_BASE) -O3 -g
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0 -DLINUX_DEBUG
CFLAGS_PROF = $(CFLAGS_BASE) -DNDEBUG -O3 -g -DPROFILE -pg -fno-inline
CFLAGS_OMP = $(CFLAGS_OPT)

FFLAGS_OPT = -O3 -g
FFLAGS_DEBUG = -g -DLINUX_DEBUG -pg -cpp -w
FFLAGS_PROF = -O3 -g
FFLAGS_OMP = $(FFLAGS_OPT)

LDFLAGS_OPT = $(LDFLAGS_BASE)
LDFLAGS_DEBUG = $(LDFLAGS_BASE)
LDFLAGS_PROF = $(LDFLAGS_BASE)
LDFLAGS_OMP = $(LDFLAGS_OPT)


ifeq ($(BUILD_MODE),OMP)
  CC += -fopenmp -DWITH_OMP
endif
