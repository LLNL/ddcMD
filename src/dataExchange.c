#include "dataExchange.h"

#include <assert.h>

#include "ddcMalloc.h"
#include "heap.h"
#include "ioUtils.h"
extern FILE *ddcfile;
/** This routine negotiates the routing for a data exchange.
 *
 *  Before calling this routine you must set nSend, destTask, and
 *  sendOffset.  These three parameters fully describe the data that
 *  will be sent.  You must also set comm.
 *
 *  This routine figures out the corresponding recv side information,
 *  nRecv, sourceTask, and recvOffset.
 *
 *  Note that the send and recv maps aren't at all generic and must be
 *  handled by some other method.
 */
void dep_negotiate(DATA_EXCHANGE_PARMS* this)
{
   const int tag = 14;
   int nTasks;
   MPI_Comm_size(this->comm, &nTasks);
   int myId;
   MPI_Comm_rank(this->comm, &myId);
   
   // How many tasks will we receive data from?
   int* msgCount = ddcMalloc(nTasks * sizeof(int));
   int* recvCnt = ddcMalloc(nTasks * sizeof(int));
   for (int ii=0; ii<nTasks; ++ii)
	{
      msgCount[ii] = 0;
      recvCnt[ii] = 1;
	}
   for (unsigned ii=0; ii<this->nSend; ++ii)
   {
      msgCount[this->destTask[ii]] += 1;
      assert(msgCount[this->destTask[ii]] == 1);
   }
   MPI_Reduce_scatter(msgCount, &this->nRecv, recvCnt, MPI_INT, MPI_SUM, this->comm);
   
   ddcFree(msgCount);
   ddcFree(recvCnt);

   // Find out which tasks we will be exchanging data with.
   int* recvBuf = (int*) ddcMalloc(2*this->nRecv*sizeof(int));
   int* sendBuf = (int*) ddcMalloc(2*this->nSend*sizeof(int));
   MPI_Request* recvReq = (MPI_Request*) ddcMalloc(this->nRecv*sizeof(MPI_Request));
   MPI_Request* sendReq = (MPI_Request*) ddcMalloc(this->nSend*sizeof(MPI_Request));
	
   for (unsigned ii=0; ii<this->nRecv; ++ii)
   {
      MPI_Irecv(recvBuf+2*ii, 2, MPI_INT, MPI_ANY_SOURCE, tag, this->comm, recvReq+ii);
   }
   for (unsigned ii=0; ii<this->nSend; ++ii)
   {
      sendBuf[2*ii] = myId;
      sendBuf[2*ii+1] = this->sendOffset[ii+1] - this->sendOffset[ii];
      int dest = this->destTask[ii];
      MPI_Isend(sendBuf+2*ii, 2, MPI_INT, dest, tag, this->comm, sendReq+ii);
   }
   MPI_Waitall(this->nSend, sendReq, MPI_STATUSES_IGNORE);
   MPI_Waitall(this->nRecv, recvReq, MPI_STATUSES_IGNORE);
	
   // set sourceTask and recvOffset arrays.
   this->sourceTask = (unsigned*) ddcMalloc(this->nRecv*sizeof(unsigned));
   this->recvOffset = (unsigned*) ddcMalloc((this->nRecv+1)*sizeof(unsigned));
   this->recvOffset[0] = 0;
   for (unsigned ii=0; ii<this->nRecv; ++ii)
   {
      this->sourceTask[ii] = recvBuf[2*ii];
      this->recvOffset[ii+1] = this->recvOffset[ii] + recvBuf[2*ii+1];
   }

   ddcFree(recvBuf);
   ddcFree(sendBuf);
   ddcFree(recvReq);
   ddcFree(sendReq);
}


/** Performs a data exchange on the sendData and recvData using the
 *  comm table information stored in this.  
 *
 *  First, the initHook function is called to perform any setup actions
 *  that may be necessary for the sendBuf and recvBuf.  Typically,
 *  initHook is a pointer to a default implementation,
 *  dep_defaultInitHook.  This default implementatation acquires memory
 *  for the sendBuf, copies data from sendData to sendBuf using sendMap
 *  and points recvBuf directly at recvData (i.e., the data will be
 *  received directly into its final storage).
 *
 *  Once the send and recv buffers are set we enter the actual
 *  communication stage: we loop through the tasks that we expect to
 *  receive data from to post receives, then through the tasks we must
 *  send to and post sends.  Two calls to MPI_Waitall ensure that all
 *  communication completes.
 *
 *  Finally, we call the freeHook to clean up.  The default
 *  implementation of the freeHook merely returns the memory used by the
 *  send buffer to the system.  However, it could be used, for example,
 *  to copy data out of the recvBuf into different storage locations if
 *  such a copy was necessary.
 */
void dep_exchange(void* sendData, void* recvData, DATA_EXCHANGE_PARMS* this)
{
   char* sendBuf = NULL;
   char* recvBuf = NULL;
   void* thunk;
   this->initHook(sendData, &sendBuf, recvData, &recvBuf, &thunk, this);
   assert(sendBuf != NULL);
   assert(recvBuf != NULL);
   
   const int tag = 15;
   MPI_Request recvReq[this->nRecv];
   MPI_Request sendReq[this->nSend];

   for (unsigned ii=0; ii<this->nRecv; ++ii)
   {
      unsigned sender = this->sourceTask[ii];
      unsigned nItems = this->recvOffset[ii+1] - this->recvOffset[ii];
      unsigned len = nItems * this->width;
      char* recvPtr = recvBuf + this->recvOffset[ii]*this->width;
      MPI_Irecv(recvPtr, len, MPI_CHAR, sender, tag, this->comm, recvReq+ii);
   }

   for (unsigned ii=0; ii<this->nSend; ++ii)
   {
      unsigned target = this->destTask[ii];
      unsigned nItems = this->sendOffset[ii+1] - this->sendOffset[ii];
      unsigned len = nItems * this->width;
      char* sendPtr = sendBuf + this->sendOffset[ii]*this->width;
      MPI_Isend(sendPtr, len, MPI_CHAR, target, tag, this->comm, sendReq+ii);
   }

   MPI_Waitall(this->nSend, sendReq, MPI_STATUS_IGNORE);
   MPI_Waitall(this->nRecv, recvReq, MPI_STATUS_IGNORE);

   this->freeHook(sendData, &sendBuf, recvData, &recvBuf, &thunk, this);
   assert(sendBuf == NULL);
   assert(recvBuf == NULL);
}




void dep_free(DATA_EXCHANGE_PARMS* this)
{
   ddcFree(this->destTask);
   ddcFree(this->sourceTask);
   ddcFree(this->sendOffset);
   ddcFree(this->recvOffset);
   ddcFree(this->sendMap);
   ddcFree(this->recvMap);
}

unsigned dep_nRemote(DATA_EXCHANGE_PARMS d)
{
   return d.recvOffset[d.nRecv];
}

typedef struct DefaultHookParms_st
{
   unsigned sendBlk;
} DEFAULT_HOOK_PARMS;

void dep_defaultInitHook(void* sendDataVoid, char** sendBuf,
                         void* recvData, char** recvBuf,
                         void** thunk, struct DataExchangeParms_st* this)
{
   DEFAULT_HOOK_PARMS* parms = ddcMalloc(sizeof(DEFAULT_HOOK_PARMS));
   *thunk = parms;
   char* sendData = sendDataVoid;
   *sendBuf = heapGet(&(parms->sendBlk));
   unsigned sendBufItems = this->sendOffset[this->nSend];
   heapEndBlock(parms->sendBlk, sendBufItems*this->width);

   for (unsigned ii=0; ii<sendBufItems; ++ii)
   {
      copyBytes((*sendBuf)+ii*this->width,
                sendData+this->sendMap[ii]*this->width,
                this->width);
   }
   
   *recvBuf = recvData;
}

void dep_defaultFreeHook(void* sendData, char** sendBuf,
                         void* recvData, char** recvBuf,
                         void** thunk, struct DataExchangeParms_st* this)
{
   DEFAULT_HOOK_PARMS* parms = *thunk;
   heapFree(parms->sendBlk);
   *sendBuf=NULL;
   *recvBuf=NULL;
   ddcFree(parms);
}


unsigned *dep_getSourceTask(DATA_EXCHANGE_PARMS* this,unsigned *nRecv)
{
	*nRecv = this->nRecv; 
	return this->sourceTask; 
}
int dep_particleIndexToSourceTask(DATA_EXCHANGE_PARMS* this,unsigned remoteParticleIndex)
{
	
	int task =-1; 
	if (this->recvMap != NULL) 
	{
	}
	for (unsigned itask =0; itask<this->nRecv;itask++)
	{
		if (remoteParticleIndex < this->recvOffset[itask+1]) 
		{
			task = this->sourceTask[itask]; 
			break ; 
		}
	}
	return task; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
