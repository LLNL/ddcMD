#ifndef DATA_EXCHANGE_H
#define DATA_EXCHANGE_H

#include <mpi.h>

/** Implements a generic data exchange framework.
 *
 *  LIMITATIONS:
 *
 *  This implementation strives for simplicity of implementation at the
 *  expense of optimization.  One consequence is that this method uses
 *  more memory for intermediate storage that is necessary for at least
 *  some cases.
 *
 *  The send and receive buffers are assumed to composed of
 *  non-overlapping sections.  This means whem multiple copies of the
 *  same data are sent to multiple receivers there will be several
 *  copies in the send buffer.
 *
 *  All sends and receives are processed together.  This framework does
 *  not exploit potential memory savings from doing only one send or
 *  recv at a time to reduce intermediate storage space.
 */


struct DataExchangeParms_st;

typedef void (*dep_exchangeHook) (void* sendData, char** sendBuf,
				  void* recvData, char** recvBuf,
				  void** thunk, struct DataExchangeParms_st* this);


typedef struct DataExchangeParms_st
{
   unsigned width;       // nBytes per data item.
   unsigned nSend;       // nTasks this task sends data to
   unsigned nRecv;       // nTasks this task recvs data from
   unsigned* destTask;   // ranks to send to
   unsigned* sourceTask; // ranks to recv from
   unsigned* sendOffset; // offsets to sendBuf
   unsigned* recvOffset; // offsets to recvBuf
   unsigned* sendMap;    // map to copy data to sendbuf
   unsigned* recvMap;    // map to copy data from recvbuf
   MPI_Comm comm;
   dep_exchangeHook initHook;
   dep_exchangeHook freeHook;
} DATA_EXCHANGE_PARMS;


void dep_negotiate(DATA_EXCHANGE_PARMS* this);
void dep_exchange(void* sendData, void* recvData, DATA_EXCHANGE_PARMS* this);
void dep_free(DATA_EXCHANGE_PARMS* this);
unsigned dep_nRemote(DATA_EXCHANGE_PARMS d);

void dep_defaultInitHook(void* sendData, char** sendBuf,
			 void* recvData, char** recvBuf,
			 void** thunk, struct DataExchangeParms_st* this);
void dep_defaultFreeHook(void* sendData, char** sendBuf,
			 void* recvData, char** recvBuf,
			 void** thunk, struct DataExchangeParms_st* this);
unsigned *dep_getSourceTask(DATA_EXCHANGE_PARMS* this,unsigned *nRecv);
int dep_particleIndexToSourceTask(DATA_EXCHANGE_PARMS* this,unsigned particleIndex);
#endif // #ifndef DATA_EXCHANGE_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
