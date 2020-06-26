#ifndef ROUTINE_H
#define ROUTINE_H
#include <mpi.h>
#include "commandLineOptions.h"

#define maxDataSize  1024

typedef struct intPair_st
{
   int first;
   int second;
} INTPAIR;


typedef struct routine_st
{
   char *name;		/* name of the routine */
   char *objclass;  
   void (*fcn)(void *, MPI_Comm);
   void *parms;
   int master; 
   int ntasks;
   int beginRank; // rank in MPI_COMM_WORLD of root rank for this routine
   int endRank; // one past the end (like a C++ iterator)
   MPI_Comm comm; 
} ROUTINE;
typedef struct routineManager_st
{
   char *name;		/* name of the system */
   char *objclass;
   char *value;
   char *type;		/* type label */
   void *parent; 
   int nRoutines;   /* number of different type of responsibilities */
   int overlap;    // should comms overlap?
   ROUTINE** routines; 
} ROUTINEMANAGER;

ROUTINE* routineManager_init(void *parent, char* name, void (* masterFcn)(void *masterParms,MPI_Comm Comm), void *masterParms);
ROUTINE *routineAllocate(int nTasks,  const char* name);
void routineDestroy(ROUTINE *routine); 
ROUTINE *routineFind(ROUTINEMANAGER *routineManager, const char* name);
MPI_Comm routineManager_getComm(ROUTINEMANAGER* routineManager, const char* name);
int routineManager_getNtasks(ROUTINEMANAGER*, const char* name); 
int routineManager_getWorldRoot(ROUTINEMANAGER *routineManager, const char* name);
INTPAIR routineManager_getWorldRankRange(ROUTINEMANAGER* routineManager, const char* name);
int *routineManager_getWorldRankList(ROUTINEMANAGER*, const char* name); 
void slaveReturn(void); 
void  callRoutine(char* name, void (*fcn)(),int dataLength, void *data);
void routineReturnSlaves(void); 
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
