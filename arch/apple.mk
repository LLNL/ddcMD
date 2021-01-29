# $Id: makefile.h 1149 2009-12-18 17:35:51Z draeger1 $ 
# Make configuration for Apple OSX

ARCHDESC = Apple OSX, mpicc/mpif90 compilers

#set this variable if an include file is used for extra dependencies
ARCHINCLUDE = 1

CPP = mpicxx -std=c++11 
CC = mpicc -std=gnu99


MACPORTS_ROOT = $(HOME)/inst/macports
#MACPORTS_ROOT = /opt/local/

# FFT library needs to be included here
CFLAGS_BASE =  -DUSE_FFTW2 -I$(MACPORTS_ROOT)/include -DWITH_MPI -DWITH_PIO -D_GNU_SOURCE -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas
LDFLAGS_BASE = -lm -lc -L$(MACPORTS_ROOT)/lib -ldfftw -lstdc++ -lquadmath

CXXFLAGS_BASE = $(CFLAGS_BASE)
CXX11FLAGS = -std=c++11 $(CFLAGS_BASE)

HAVE_GSL = 0
ifeq ($(HAVE_GSL),0)
CFLAGS_BASE  += -DHAVE_GSL
LDFLAGS_BASE += -lgsl -lgslcblas
endif

HAVE_QHULL = 0
ifeq ($(HAVE_QHULL),1)
  # MacPorts version of qhull (2010.1) is incompatible with 
  # current qhullInterface.h 
  #QHULLLIB        = $(MACPORTS_ROOT)/lib/libqhull.a
  #QHULLINCLUDEDIR = $(MACPORTS_ROOT)/include/qhull
  LIBQHULL        = ../base/qhull/libqhull_apple.a
  QHULLINCLUDEDIR = ../base/qhull/src
  CFLAGS_BASE += -DUSE_QHULL -I$(QHULLINCLUDEDIR)
endif

CFLAGS_OPT = $(CFLAGS_BASE) -O3
CFLAGS_DEBUG = $(CFLAGS_BASE) -g -O0 -DLINUX_DEBUG
CFLAGS_PROF = $(CFLAGS_BASE) -O3 -DPROFILE -pg

#erik: debug and prof seem entangled, but this is what was in makefile.h
FFLAGS_OPT = -O3
FFLAGS_DEBUG = -g -DLINUX_DEBUG -pg -cpp -w
FFLAGS_PROF = -O3

LDFLAGS_OPT = $(LDFLAGS_BASE)
LDFLAGS_DEBUG = $(LDFLAGS_BASE)
LDFLAGS_PROF = $(LDFLAGS_BASE)
