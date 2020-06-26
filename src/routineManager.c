#include "routineManager.h"
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "object.h"
#include "ddcMalloc.h"
#include "external.h"
#include "commandLineOptions.h"
#include "mpiUtils.h"
#include "utilities.h"
#include "commInfo.h"
#include "ioUtils.h"
#include "ptiming.h"

#define MAX(A, B) ((A) > (B) ? (A) : (B))

extern MPI_Comm COMM_LOCAL;

typedef struct rendezvous_msg_str
{
   int action;
   void (*fcn)();
   int dataLength;
   char  data[maxDataSize];
} RENDEZVOUSMSG;

typedef struct descriptor_st
{
   char *name;
   char *objclass;
} DESCRIPTOR;  


static int _nextAvailableRank =0; 
static ROUTINEMANAGER *_routineManager = NULL; 


static ROUTINE* createRoutine(int ntasks, int overlap);

/**
 *  When overlap is specified a task will belong to more than one
 *  routine.  In that case this routine needs to choose which routine to
 *  return.  We call the routine that will be returned the primary
 *  routine and define it as routine with the largest number of tasks
 *  associated with it.  If multiple routines have the same number of
 *  tasks then we choose the member of that set that is listed first in
 *  the routines keyword.
 *
 *  This function also sets COMM_LOCAL. 
 */
ROUTINE* routineManager_init(void *parent, char* name, void (* masterFcn)(void *parms, MPI_Comm Comm), void *masterParms)
{
   timestampBarrier("Start routineManager", MPI_COMM_WORLD);
   assert(_routineManager == NULL);
   _routineManager = (ROUTINEMANAGER *) object_initialize(name, "ROUTINEMANAGER", sizeof(ROUTINEMANAGER));
   
   //DESCRIPTOR *routineDescriptor;  
   int layouts[16]; 
   int nLayouts =  object_get((OBJECT*)_routineManager, "layouts",  layouts,   INT,16,"-1");   
   //int nRoutines = object_getv((OBJECT*)_routineManager, "routines",  (void*)&routineDescriptor,   STRING,IGNORE_IF_NOT_FOUND);   
   object_get((OBJECT*)_routineManager, "overlap", &_routineManager->overlap, INT, 1, "0");

   //nRoutines/=2;
   int nRoutines=nLayouts; 
   if (nLayouts != nRoutines)
   {
      printf("ERROR:  Error in ROUTINEMANAGER object.\n"
				 "        Number of layouts does not match number of routines.\n");
      abortAll(-1);
   }
   
   _routineManager->nRoutines = nRoutines; 
   _routineManager->routines = ddcMalloc(nRoutines*sizeof(ROUTINE*)); 
   if (layouts[0] == -1 ) // runs on all available tasks
      layouts[0]=getSize(-1);

   if (! _routineManager->overlap)
   {
      int layoutSum = 0;
      for (int ii=0; ii<nRoutines; ++ii)
			layoutSum += layouts[ii];
      if (layoutSum != getSize(-1) && getRank(0) == 0)
      {
			printf("ERROR:  Error in ROUTINEMANAGER object\n"
					 "        Sum of layouts does not match number of tasks in job.\n");
			abortAll(-1);
      }
   }

   for (int i=0; i<nRoutines; i++)
   {
      WAIT(-1); 
      ROUTINE* routine = _routineManager->routines[i] = createRoutine(layouts[i], _routineManager->overlap); 
      if (i==0) 
      { 
         routine->master=1;
         routine->fcn=masterFcn;
         routine->parms = masterParms; 
      }
      if (routine->comm == MPI_COMM_NULL) continue;
      int myRank;
      MPI_Comm_rank(routine->comm, &myRank);
      if (myRank == 0)
      {
			int isRect = isRectangularComm(routine->comm);
			if (isRect == 0)
				printf("  Communicator for routine # %d is NOT rectangular\n", i);
			else
				printf("  Communicator for routine # %d is rectangular\n", i);
      }
   }

   // Find primary routine for task.
   ROUTINE* primaryRoutine = NULL; 
   for (int ii=0; ii<nRoutines; ++ii)
   {
      if (_routineManager->routines[ii]->comm == MPI_COMM_NULL) continue;
      if (primaryRoutine != NULL  &&
			 primaryRoutine->ntasks >= _routineManager->routines[ii]->ntasks)
			continue;
      primaryRoutine = _routineManager->routines[ii];
   }

   if (primaryRoutine == NULL)
   {
      printf("ERROR:  Error in ROUTINEMANAGER object\n"
	     "        Task %d not assigned to any routine.\n", getRank(-1));
      abortAll(-1);
   }
   COMM_LOCAL = primaryRoutine->comm; 
   timestampBarrier("Finish routineCreate", MPI_COMM_WORLD);
   ROUTINE** routines = _routineManager->routines; 
   routines[0]->name =strdup("routineMaster"); 
   
   return primaryRoutine; 
}

void routineDestroy(ROUTINE *routine) 
{
   if (routine->master) routineReturnSlaves();
   COMM_LOCAL = MPI_COMM_WORLD; 
}

ROUTINE *routineFind(ROUTINEMANAGER *routineManager, const char* name)
{
   if (routineManager==NULL) routineManager = _routineManager;
   for (int i=0;i<routineManager->nRoutines;i++)
   {
      ROUTINE *routine=routineManager->routines[i]; 
	  if (routine->name == NULL) continue; 
      if (strcmp(routine->name,name)==0) return routine;
   }
   return NULL;
}
MPI_Comm routineManager_getComm(ROUTINEMANAGER *routineManager, const char* name)
{
	ROUTINE *routine = routineFind(routineManager,name); 
    if (routine == NULL) return MPI_COMM_NULL; 
    return routine->comm; 
}

int routineManager_getNtasks(ROUTINEMANAGER *routineManager, const char* name)
{
	ROUTINE *routine = routineFind(routineManager,name); 
   if (routine ==NULL) return 0;
    return routine->ntasks;
}

int routineManager_getWorldRoot(ROUTINEMANAGER *routineManager, const char* name)
{
	ROUTINE *routine = routineFind(routineManager, name); 
   if (routine == NULL) return -1;
   return routine->beginRank;
}

INTPAIR routineManager_getWorldRankRange(ROUTINEMANAGER* routineManager, const char* name)
{
   INTPAIR tmp;
   tmp.first = 0;
   tmp.second = 0;
   ROUTINE *routine = routineFind(routineManager, name); 
   if (routine != NULL)
   {
      tmp.first = routine->beginRank;
      tmp.second = routine->endRank;
   }
   return tmp;
}

void  callRoutine(char* name, void (*fcn)(),int dataLength, void *data)
{

   assert(dataLength<=maxDataSize); 
   RENDEZVOUSMSG msg; 
   int masterLeaderRank, slaveLeaderRank;
   int rendezvous_tag = 99;
   int commLocalRank;
   msg.action=1; 
   msg.fcn=fcn;
   msg.dataLength=dataLength;
   copyBytes(msg.data,data,dataLength); 
   if (fcn == slaveReturn ) msg.action=0;

   MPI_Comm_rank( COMM_LOCAL, &commLocalRank );
   if (commLocalRank !=0) return; 

   masterLeaderRank = _routineManager->routines[0]->beginRank ; // Get the rank of the leader of the master group in MPI_COMM_WORLD
   if (name != NULL) 
   {
      slaveLeaderRank = routineManager_getWorldRoot(NULL,name); // Get the rank of the leader of the slave group in MPI_COMM_WORLD
      if (slaveLeaderRank == masterLeaderRank)  return; 

      if ( slaveLeaderRank >=0 ) MPI_Send( (void*)(&msg), sizeof(msg), MPI_BYTE, slaveLeaderRank, rendezvous_tag, MPI_COMM_WORLD); // Send a message to the slave leader (if there is one)
   }
   else 
   { 
      for (int i=1;i<_routineManager->nRoutines;i++)
      { 
         ROUTINE *routine=_routineManager->routines[i]; 
         slaveLeaderRank = routine->beginRank; 
         if (slaveLeaderRank == masterLeaderRank)  continue; 
         if ( slaveLeaderRank >=0 ) MPI_Send( (void*)(&msg), sizeof(msg), MPI_BYTE, slaveLeaderRank, rendezvous_tag, MPI_COMM_WORLD); // Send a message to the slave leader (if there is one)
      }
   } 
}
void  routineReturnSlaves(void)
{

   RENDEZVOUSMSG msg; 
   int masterLeaderRank, slaveLeaderRank;
   int rendezvous_tag = 99;
   int commLocalRank;
   msg.action=0; 
   msg.fcn=NULL;
   msg.dataLength=0;

   MPI_Comm_rank( COMM_LOCAL, &commLocalRank );

   masterLeaderRank = _routineManager->routines[0]->beginRank ; // Get the rank of the leader of the master group in MPI_COMM_WORLD

   for (int i=1;i<_routineManager->nRoutines;i++) 
   { 
      slaveLeaderRank = _routineManager->routines[i]->beginRank ; // Get the rank of the leader of the master group in MPI_COMM_WORLD
      if (slaveLeaderRank == masterLeaderRank)  continue; 
      if ( commLocalRank == 0  && slaveLeaderRank >=0 ) MPI_Send( (void*)(&msg), sizeof(msg), MPI_BYTE, slaveLeaderRank, rendezvous_tag, MPI_COMM_WORLD); // Send a message to the slave leader (if there is one)
   }


}
ROUTINE *routineMaster(void)
{
   return _routineManager->routines[0]; 
}
ROUTINE *routineCurrent(void)
{
   for (int i=0;i<_routineManager->nRoutines;i++) 
   { 
      ROUTINE *routine=_routineManager->routines[i]; 
      if ( routine->comm == COMM_LOCAL) return routine; 
   }
   return NULL; 
}
void routineName(char  *data)
{
   int routineIndex = *(int *) data; 
   char *name = data+4; 
   ROUTINE *routine = _routineManager->routines[routineIndex]; 
   routine->name = strdup(name); 
}
ROUTINE *routineAllocate(int nTasks, const char* name)
{
   ROUTINE *routineMaster=routineCurrent(); 
   if (!routineMaster->master) return NULL; 
   int masterLeaderRank, slaveLeaderRank;
   int rendezvous_tag = 99;
   int commLocalRank;

   MPI_Comm_rank( COMM_LOCAL, &commLocalRank );

   masterLeaderRank = routineMaster->beginRank ; // Get the rank of the leader of the master group in MPI_COMM_WORLD

   ROUTINE *routine=NULL; 
   int routineIndex=-1; 
   for (routineIndex=1;routineIndex<_routineManager->nRoutines;routineIndex++) 
   { 
      routine=_routineManager->routines[routineIndex]; 
      if ( routine->name == NULL && ((routine->ntasks==nTasks)|(_routineManager->overlap !=0)) ) break; 
   }
   if (routine == NULL) return NULL; 
   slaveLeaderRank = routine->beginRank ; // Get the rank of the leader of the master group in MPI_COMM_WORLD
   routine->name = strdup(name); 
   if (slaveLeaderRank == masterLeaderRank)  return routine; 
   if ( commLocalRank == 0  && slaveLeaderRank >=0 ) 
   {
      RENDEZVOUSMSG msg; 
      msg.action=1; 
      msg.fcn=routineName ;
      msg.dataLength=5+strlen(name);
      copyBytes(msg.data,(char *)&routineIndex,4); 
      copyBytes(msg.data+4,name,msg.dataLength-4); 
      MPI_Send( (void*)(&msg), sizeof(msg), MPI_BYTE, slaveLeaderRank, rendezvous_tag, MPI_COMM_WORLD); // Send a message to the slave leader (if there is one)
   }
   return routine; 

}

int listenRoutine(RENDEZVOUSMSG *msg)
{

   int masterLeaderRank;
   int rendezvous_tag = 99;
   int commLocalRank;
   MPI_Status status;

   profile(RENDEZVOUS_T, START);
   MPI_Comm_rank( COMM_LOCAL, &commLocalRank );

   masterLeaderRank = _routineManager->routines[0]->beginRank ; // Get the rank of the leader of the master group in MPI_COMM_WORLD
   if ( commLocalRank == 0 ) MPI_Recv( msg, sizeof(*msg), MPI_BYTE, masterLeaderRank, rendezvous_tag, MPI_COMM_WORLD, &status ); // Wait for a message from master leader
   MPI_Bcast( msg, sizeof(*msg), MPI_BYTE, 0, COMM_LOCAL ); //share message with group members

   profile(RENDEZVOUS_T, END);
   return msg->action;
}

void slave(void *parms, MPI_Comm Comm)
{

   RENDEZVOUSMSG msg; 
   while (listenRoutine(&msg))
   {
      msg.fcn(msg.data);    
   }
   /*
   char *name; 
   if (routine->name != NULL) name = routine->name; else name ="unallocated slave"; 
   char string[10+strlen(name)]; 
   sprintf(string,"Finished %s",name); 
   */

   timestampBarrier("Finished slave",Comm);
}
void  slaveReturn (void) { }


/** When overlap is non-zero the comm for this task will start with
 * world rank zero regardless of the value of _nextAvailableRank.  This
 * will cause the comms from different routines to overlap.  */
ROUTINE* createRoutine(int nTasks, int overlap)
{   
   int worldRank, worldSize; 
   MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
   MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
   ROUTINE *routine = ddcMalloc(sizeof(ROUTINE)); 
   routine->master=0; 

   int beginRank = _nextAvailableRank;
   if (overlap != 0)
      beginRank = 0;
   int endRank = beginRank + nTasks;

   assert(endRank <= worldSize);
   
   _nextAvailableRank = MAX(_nextAvailableRank, endRank);

   MPI_Comm comm=MPI_COMM_NULL; 
   if (nTasks == worldSize)
   {
      comm = MPI_COMM_WORLD;
   }
   else
   {
      int color = MPI_UNDEFINED;
      if (worldRank >= beginRank && worldRank < endRank)
         color = 1;
      
      timestampBarrier("Start comm_split in createRoutine",MPI_COMM_WORLD);
      MPI_Comm_split(MPI_COMM_WORLD, color, worldRank, &comm);
      timestampBarrier("Finished comm_split in createRoutine",MPI_COMM_WORLD);
      
      if (color == MPI_UNDEFINED) 
         comm = MPI_COMM_NULL; // This is what MPI_Comm_split should return anyway.
   }
   
   routine->ntasks = nTasks; 
   routine->beginRank = beginRank;
   routine->endRank = endRank;
   routine->comm = comm;
   routine->fcn = slave; 
   routine->parms = NULL; 
   routine->name =NULL; 
   return routine; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
