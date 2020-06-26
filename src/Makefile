#$Id: $  
#########################################################
# PARALLELIZED MAKEFILE
# This is a minially invasive change of the Makefile. 
# All the magic happens HERE. No other files have been 
# modified (see the commit).
#########################################################

SHELL = /bin/bash

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o .a

### .SUFFIXES: .a

.PHONY: DEFAULT opt debug profile makedirs clean distclean depend install

# Get name of architecture-specific include file from hostname
USE_GPU=0
include ../arch/Makefile.arch

# If defined, use architecture file set with ARCH variable
ARCH ?= $(ARCHGUESS)
ARCHFILE = ../arch/$(ARCH).mk
include $(ARCHFILE)

ifdef DEVLPMT 
include devlpmt/devlpmt.mk
CFLAGS_BASE += -DDEVLPMT 
LIST := $(shell ./listfile *devlpmt.* )
endif
#GITVERSION := $(shell git rev-parse HEAD )
GITVERSION := $(shell git log --pretty=format:"%h [%cd]" -1 )

CFLAGS_TEST = $(CFLAGS_BASE) -DTESTEXE=1
CFLAGS_TESTD =  $(CFLAGS_BASE) -DTESTEXE=1 -g -O0 -DLINUX_DEBUG

BINDIR  =  ../bin
LIBDDC = $(LIBDIR)/lib$(EXENAME).a
LIBDDCCPP = $(LIBDIR)/lib$(EXENAME)CPP.a
ifeq ($(USE_GPU),1)
LIBDDCGPU = $(LIBDIR)/lib$(EXENAME)GPU.a
endif

DEFAULT: opt

opt:
	@echo "Building ddcMD with $@ options, using $(ARCHFILE) ($(ARCHDESC))"
	@echo "Building with GPU support: "$(USE_GPU)
	@echo "..."
	if [ ! -z $(MODULES) ]; then module load $(MODULES); fi
	$(MAKE) makedirs $(BINDIR)/ddcMD-$(ARCH) CFLAGS="$(CFLAGS_OPT)" \
	  NVCCFLAGS="$(NVCCFLAGS_OPT)" \
	  FFLAGS="$(FFLAGS_OPT)" LDFLAGS="$(LDFLAGS_OPT)" \
	  LIBDIR=../lib/$(ARCH)-opt EXENAME=ddcMD-$(ARCH) DEPBASE=$(ARCH)_opt

debug:
	@echo "Building ddcMD with $@ options, using $(ARCHFILE) ($(ARCHDESC))"	
	@echo "Building with GPU support: "$(USE_GPU)
	@echo "..."
	if [ ! -z $(MODULES) ]; then module load $(MODULES); fi
	@$(MAKE) --no-print-directory makedirs $(BINDIR)/ddcMD-$(ARCH)-debug CFLAGS="$(CFLAGS_DEBUG)" \
	  NVCCFLAGS="$(NVCCFLAGS_DEBUG)" \
	  FFLAGS="$(FFLAGS_DEBUG)" LDFLAGS="$(LDFLAGS_DEBUG)" \
	  LIBDIR=../lib/$(ARCH)-debug EXENAME=ddcMD-$(ARCH)-debug DEPBASE=$(ARCH)_debug
profile:
	@echo "Building ddcMD with $@ options, using $(ARCHFILE) ($(ARCHDESC))"
	@echo "Building with GPU support: "$(USE_GPU)
	@echo "..."
	if [ ! -z $(MODULES) ]; then module load $(MODULES); fi
	@$(MAKE) --no-print-directory makedirs $(BINDIR)/ddcMD-$(ARCH)-prof CFLAGS="$(CFLAGS_PROF)" \
	  NVCCFLAGS="$(NVCCFLAGS_PROF)" \
	  FFLAGS="$(FFLAGS_PROF)" LDFLAGS="$(LDFLAGS_PROF)" \
	  LIBDIR=../lib/$(ARCH)-prof EXENAME=ddcMD-$(ARCH)-prof DEPBASE=$(ARCH)_prof
omp:
	module load $(MODULES)
	@echo "Building ddcMD with $@ options, using $(ARCHFILE) ($(ARCHDESC))"
	@echo "Building with GPU support: "$(USE_GPU)
	@echo "..."
	if [ ! -z $(MODULES) ]; then module load $(MODULES); fi
	  CFLAGS="$(CFLAGS_OMP)" FFLAGS="$(FFLAGS_OMP)" LDFLAGS="$(LDFLAGS_OMP)" \
	  LIBDIR=../lib/$(ARCH)-omp EXENAME=ddcMD-$(ARCH)-omp DEPBASE=$(ARCH)_omp
install: 
	cp $(BINDIR)/ddcMD-$(ARCH) $(INSTALL_DIR) 

notUsedSRC    = forcetest.c unionGroup.c 

miscSRC      = eq.c threadScheduler.c functions.c pair_function.c saveState.c parityHandler.c bglParity.c pinfo.c 

SRC          = ddcMD.c commons.c routineManager.c simulate.c masterFactory.c masters.c system.c  molecule.c molecularPressure.c  species.c  box.c  boxPrescriptiveTime.c state.c \
               collection.c collection_write.c collection_read.c  energyInfo.c integrationTest.c testingUtils.c forcetest.c accelerator.c


recbislibSRC =   alltoall_sparse.c  balance.c       equalize.c  genptr.c   kstat.c     rectimer.c      treecomm.c \
                 approx_median.c    balance_hist.c  exchange.c  integer.c  parkstat.c  redistribute.c

recbisddcmdSRC = state_redist.c back_communicate.c bisectionLoadBalance.c recnlist.c  bisectionCalc.c

recbisSRC      = $(recbislibSRC) $(recbisddcmdSRC)

ddcSRC       = ddcFindPairs.c ddcUpdate.c ddcUpdateAll.c sortFunctions.c ddcenergy.c \
               domain.c md2ddc.c ddc.c ddcRule.c ddcRuleMolecule.c ddcAssignment.c ddcSendRecv.c  ddc_writePXYZ.c readPXYZ.c particle.c 

loadBalanceSRC = loadBalance.c zRampLoadBalance.c  

bigFFTSRC    = ProcessGrid2D.c 

mathSRC      = crc32.c tfunction.c spline.c erfErfc.c  gsl.c  sph.c \
               ErfcRootXOverRootX.c ExpRootXOverRootX.c  primes.c lcg64.c random.c solve.c tridiagLU.c   $(bigFFTSRC)

utilSRC       = recordDecorator.c gidShuffle.c gid.c commInfo.c mpiAlgorithm.c mpiTypes.c\
                ompUtils.c dataExchange.c  dataExchangeUtils.c getRemoteData.c \
                qIntrinsics.c makeImageParticle.c

memorySRC    = expandbuffer.c resize.c 

simUtilSRC   = GridAssignmentObject.c error.c heap.c ioUtils.c match.c object.c pioFixedRecordHelper.c pioVariableRecordHelper.c three_algebra.c utilities.c\
               ddcMalloc.c hardwareInfo.c intQueue.c lessThan.c mpiUtils.c pio.c pioHelper.c units.c tagServer.c

ioSRC        = io.c sCatPrintf.c objectSetup.c check_line.c readCmds.c commandLineOptions.c  printinfo.c graph.c \
               profile.c  hpmWrapper.c format.c

groupSRC     = group.c berendsen.c doubleMirror.c extforce.c fixedVelocity.c free.c frozen.c langevin.c\
               piston.c quench.c shear.c \
               shwall.c  shock.c  unionGroup.c 

pairSRC      = neighbor.c pairlist.c geom.c preduce.c auxNeighbor.c \
               pairIterator.c pairFinder.c pairFinderGeom.c pairFinderModified.c rtuPairFinder.c nearPair.c

integratorSRC = integrator.c nextfile.c  nptglf.c nglf.c nglfNew.c nglfNK.c nglfrattle.c nglfconstraint.c nglfError.c nglfTest.c nveglf.c  nveglf_simple.c  


analysisSRC  = analysis.c ackland_jones.c  centrosym.c coarsegrain.c dataSubset.c dsf.c\
               forceAverage.c kineticEnergyDistn.c paircorrelation.c pairAnalysis.c  \
               quaternion.c ssf.c stressWrite.c subsetWrite.c vcmWrite.c velocityAutocorrelation.c zdensity.c\
					cholAnalysis.c 

transformSRC = transform.c alchemyTransform.c boxTransform.c gidShuffleTransform.c customTransform.c projectileTransform.c  \
                    selectSubset.c addVelocity.c linearisotropicv.c thermalize.c thermalizeTransform.c replicate.c \
                 transectMorph.c impactTransform.c

potentialSRC  = potential.c energy.c orderSH.c pair.c pairEnergy.c pairKernel.c reflect.c \
               table_function.c zeroPotential.c

pppmSRC      = 
               

mgptSRC      = 
                


eamSRC = eam.c eam_exp.c eam_fs.c eam_sc.c eam_at.c eam_rational.c eam_tabular.c 

plasmaSRC    = 

cppSRC       = binProcess.cpp
 
gpuSRC       = minMaxGPU.cu pairProcessGPU.cu radix_sort.cu gpuMemUtils.cu prefixScan.cu nlistGPU.cu bondedGPU.cu \
               nglfGPU.cu molecularPressureGPU.cu \
               restraintGPU.cu nglfconstraintGPU.cu kineticGPU.cu bioMartini.cu runtimeKernel.cu

bioSRC       = bioCharmmTop.c bioCharmmPar.c bioString.c bioCharmmParms.c bioCharmm.c\
               bioCharmmCovalent.c bioCharmmCovalentEnergies.c  bioCharmmCovalentEnergiesNormal.c\
               bioCharmmCovalentEnergiesSorted.c bioCharmmCovalentEnergiesSortedWeighted.c \
               bioTransform.c bioCharmmRule.c bioMartini.c bioMartiniRule.c bioMMFF.c restraint.c

cdbSRC       = 

unitTestSRC =      
ifeq ($(USE_GPU),1)
unitTestSRC+=
else
pairSRC += constructList.c 
endif

//integrTestSRC = integrationTest.c


fmmSRC       = 


ALLSRC = $(SRC) $(simUtilSRC) $(ioSRC) $(utilSRC) $(miscSRC) $(recbisSRC) \
         $(ddcSRC) $(loadBalanceSRC) $(memorySRC) $(mathSRC) $(groupSRC) $(integratorSRC) $(analysisSRC)  $(transformSRC) \
         $(pairSRC) $(potentialSRC) $(pppmSRC) $(mgptSRC) $(eamSRC) $(plasmaSRC) \
         $(ARCH_SRC) $(DEVLPMT_SRC) $(bioSRC) $(fmmSRC) $(cdbSRC) $(unitTestSRC) $(integrTestSRC) nullRoutines.c 

releaseFiles: 
	ls -l $(ALLSRC) $(gpuSRC) $(cppSRC) constructList.c version.c | fgrep -v '>'  | awk '{print $$9}' | sort 
nullRoutines:
	make |& fgrep  ': undefined ref' | awk '{printf("void %s(){}\n", substr($$5,2,length($$5)-2))}'  >> nullRoutines.c 

RELEASE:
	./release.csh  "$(RELEASESRC)"
SRCEXTRA =   
ALLSRC += $(SRCEXTRA)




# Create all the compilation rules for each source type
# $(LIBDIR)/%.f.o: %.f
# 	+$(FC) $(FFLAGS) -DFMM_DOUBLE -DFMM_I32 -Ifmm/include -I./ -I../ -Lfmm/ -c -o $@ $<

$(LIBDIR)/%.o: %.c
	+$(CC) $(CFLAGS) -c -o $@ $<

$(LIBDIR)/%.cc.o: %.cc
	+$(CPP) $(CPPFLAGS_BASE) -c -o $@ $<

$(LIBDIR)/%.cpp.o: %.cpp
	+$(CPP) $(CPPFLAGS) -c -o $@ $<

ifeq ($(USE_GPU),1)
$(LIBDIR)/%.cu.o: %.cu
	+$(NVCC) $(NVCCFLAGS)  -c -o $@ $<
endif

# Create the list of compile objects to be archived
LIBDDC_OBJ    = $(ALLSRC:%.c=$(LIBDIR)/%.o) $(CPPSRC:%.cc=$(LIBDIR)/%.cc.o) 
# $(FSRC:%.f=$(LIBDIR)/%.f.o)
LIBDDCCPP_OBJ = $(cppSRC:%.cpp=$(LIBDIR)/%.cpp.o)
ifeq ($(USE_GPU),1)
LIBDDCGPU_OBJ = $(gpuSRC:%.cu=$(LIBDIR)/%.cu.o)
endif

# Appends all the lists coming from different source types
ALLOBJS = $(LIBDDC_OBJ) $(LIBDDCCPP_OBJ) 
ifeq ($(USE_GPU),1)
ALLOBJS += $(LIBDDCGPU_OBJ)
endif

# Archives the source types into separates static libraries
comp_OBJ: $(ALLOBJS) 
	+ar rv $(LIBDDC) $(LIBDDC_OBJ)

compCPP_OBJ: $(ALLOBJS)
	+ar rv $(LIBDDCCPP) $(LIBDDCCPP_OBJ)

ifeq ($(USE_GPU),1)
compGPU_OBJ: $(ALLOBJS)
	+ar rv $(LIBDDCGPU) $(LIBDDCGPU_OBJ)
endif

# Compiles libqhull
ifdef LIBQHULL
$(LIBQHULL):
	cd ../base; ./doinstall
endif

# Creates a list of all static libraries
ALL_ARCHIVES = comp_OBJ compCPP_OBJ
ifeq ($(USE_GPU),1)
ALL_ARCHIVES += compGPU_OBJ
endif

# Compiles ddcMD. Be careful when changing this line for the linking order of the libraries
#$(BINDIR)/$(EXENAME): version.c $(LIBQHULL) $(ALL_ARCHIVES)
#	$(CC) $(CFLAGS) -o $@ version.c $(LIBDDC) $(LIBDDCCPP) $(LIBDDCGPU) $(LIBQHULL) $(LDFLAGS) \
#-DDEVLPMT_SRC="$(LIST)" \
#-DGITVERSION="$(GITVERSION)" -DCFLAGS="$(CFLAGS)" \
#-DLDFLAGS="$(LDFLAGS)" -DPATH="$(PWD)" -DTARGET=$(EXENAME) -DARCH=$(ARCH)
$(BINDIR)/$(EXENAME): version.c $(LIBQHULL) $(ALLOBJS)
	$(CC) $(CFLAGS) -o $@ version.c $(ALLOBJS) $(LIBQHULL) $(LDFLAGS) \
-DDEVLPMT_SRC="$(LIST)" \
-DGITVERSION="$(GITVERSION)" -DCFLAGS="$(CFLAGS)" \
-DLDFLAGS="$(LDFLAGS)" -DPATH="$(PWD)" -DTARGET=$(EXENAME) -DARCH=$(ARCH) 
#
# Adding a check that the executable file has been successfully generated.
	@if [ -a $(BINDIR)/$(EXENAME) ]; then  echo ; echo; echo COMPILATION SUCCESSFUL ; fi
	ls $(PWD)/$(BINDIR)

makedirs: 
	@if [ ! -d ../lib ]; then mkdir -p ../lib ;  fi
	@if [ ! -d ../lib/$(ARCH)/fmm ]; then mkdir -p ../lib/$(ARCH)/fmm; fi
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR) ;  fi
	@if [ ! -d $(LIBDIR)/fmm ]; then mkdir -p $(LIBDIR)/fmm; fi
	@if [ ! -d ../bin ]; then mkdir -p ../bin ; fi

clean :
	rm -fr ../lib/$(ARCH)-*/*.a ../lib/$(ARCH)-*/*.o ../base/qhull/libqhull_$(ARCH).a
	rm -fr $(BINDIR)/ddcMD*$(ARCH)*
	rm -f .depend_*
	rm -f *.lst

distclean : clean
	rm -f $(BINDIR)/*
	rm -rf ../lib/* ../base/libqhull_*.a ../base/qhull/libqhull_*.a
	rm -f $(LIBQHULL)

AUX = Makefile Makefile.arch ../arch version.c $(wildcard *.h)  primeascii2bin.c  FORCEnoif7.h 

#version:  
#	$(CC) version.c -DCFLAGS="\"${CFLAGS}\"" -DLDFLAGS="\"${LDFLAGS}\"" -DPATH=${PWD} -o version

.depend_$(DEPBASE): $(ALLSRC) $(AUX)
	touch .depend_$(DEPBASE)
	$(MAKE) --no-print-directory depend LIBDIR=$(LIBDIR) DEPBASE=$(DEPBASE)

depend:
	echo "Rebuilding dependencies..."
	@makedepend -f .depend_$(DEPBASE) -Y. -p'$$(LIBDDC)($(LIBDIR)/' -o '.o)' --$(CFLAGS)-- $(ALLSRC) 2> /dev/null

# this suppresses deletion of intermediate .o files, but screws up dependencies
#ALLOBJS = $(ALLSRC:%.c=$(LIBDIR)/%.o) 
#ALLOBJS += $(CXXSRC:%.cc=$(LIBDIR)/%.o)
#.SECONDARY:  $(ALLOBJS)

-include ../arch/$(ARCH).include

ifdef DEPBASE
  -include .depend_$(DEPBASE)
endif
mallinfoBug: mallinfoBug.c 
	$(CC) $(CFLAGS) $(LDFLAGS_OPT) mallinfoBug.c -o mallinfoBug
