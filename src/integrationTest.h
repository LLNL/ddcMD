#ifndef INTEGRATIONTEST_H
#define INTEGRATIONTEST_H
#include "system.h" 
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "object.h"
#include "simulate.h"
#include "species.h"
#include "ddc.h"
#include "box.h"
#include "collection.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include "energyInfo.h"
#include "printinfo.h"

//
typedef struct integrationTest_st{
   char *name;		/*!< name of the integrationTest object */
   char *objclass;
   char *value;
   char *type;		/* type label */
   void  *parent; 
   MPI_Comm comm;
   
   //! number of pairs of potentials in list.  we test each potential within a pair outputs same energies. 
   //! thus each pair contains two implementations of the same potential
   int nPotentialPotentialTests;
   char **PotentialPotentialTests; /*!< list of pairs of potentials*/
 
   //! like the list of potential potential pairs, except each pair is a tuple of a potential
   //! and a data file containing correct energies for that potential
   int nPotentialDataTests;
   char **PotentialDataTests; /*!< list of pairs of potential data file pairs*/
 
   int dump_positions;    
} INTEGRATIONTEST;
INTEGRATIONTEST* integrationTest_init(void *parent, char *name, MPI_Comm comm);
int runTests(); 
int compareEnergies(ETYPE *energyInfo1, ETYPE * energyInfo2, char *potname1, char *potname2, UNITS *units,int nglobal, double tolerance );
int compareForces(STATE *s1, STATE *s2, char *potname1, char *potname2, UNITS *units,int nglobal, double tolerance );
#endif //INTEGRATIONTEST
