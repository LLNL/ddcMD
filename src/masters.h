#ifndef MASTERS_H
#define MASTERS_H
#include <mpi.h>
#include "commandLineOptions.h"
#include "simulate.h"


typedef struct commonParms_st
{
   char *simulateName; 
   char *objectFile; 
   char *restartFile; 
} COMMONPARMS;
typedef struct master_st
{
   void (*fcn)(void *parms, MPI_Comm SimulateComm);
   void *parms; 
} MASTER;
typedef struct simulateMasterParms_st
{
      COMMONPARMS common;
      time_t stopTime; 
} SIMULATEMASTERPARMS; 

typedef struct thermalizeMasterParms_st
{
      COMMONPARMS common;
      double thermalTemp; 
} THERMALIZEMASTERPARMS; 

typedef struct analysisMasterParms_st
{
      COMMONPARMS common;
      int energyCalc; 
} ANALYSISMASTERPARMS; 
typedef struct defaultMasterParms_st
{
      COMMONPARMS common;
} DEFAULTMASTERPARMS; 
MASTER masterFactory(COMMAND_LINE_OPTIONS  opt); 
void transformMaster(void *parms, MPI_Comm SimulateComm);
void eightFoldMaster(void *parms, MPI_Comm SimulateComm);
void thermalizeMaster(void *parms, MPI_Comm SimulateComm);
void analysisMaster(void *parms, MPI_Comm SimulateComm);
void readWriteMaster(void *parms, MPI_Comm SimulateComm);
void testPressureMaster(void *parms, MPI_Comm SimulateComm);
void testForceMaster(void *parms, MPI_Comm SimulateComm);
void simulateMaster(void *parms, MPI_Comm SimulateComm);
void firstEnergyCall(SIMULATE *simulate);
void integrationTestMaster(void *parms, MPI_Comm SimulateComm);
void unitTestMaster(void *parms, MPI_Comm SimulateComm);
void adjustBox(SIMULATE *simulate);
#endif
