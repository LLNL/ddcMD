#include "getRemoteData.h"
#include <mpi.h>
#include <assert.h>

#include "ddc.h"

/** Returns 1 if a is sorted in ascending order, 0 otherwise */
static int sorted(const unsigned* a, unsigned n)
{
   for (unsigned ii=1; ii<n; ++ii)
      if (a[ii] < a[ii-1])
	 return 0;
   return 1;
}

/**
 *  This function allows tasks to request remote data for specific
 *  particles.  It exists to serve the case where remote data is needed
 *  rarely and/or for only a few particles.  Using this function to get
 *  the data avoids the need to register remote data fields and drag the
 *  data around for every particle at every time step.
 *
 *  This function is not designed for situations where data is requested
 *  for a large number of particles.  There is a hard-coded limit for
 *  the largest number of requests that any remote task can receive of
 *  maxRequests = 512.
 *
 *  Arguments:
 *
 *  iLocal An array containing the local indices (in the STATE) of the
 *     particles for which we wish to obtain remote information.
 *     Obviously each index should be larger than nLocal since
 *     there is no sense in asking for data that we already have
 *     locally.  The indices must be sorted in ascending order.
 *
 *  iLocalSize The number of particles for which data is requested.
 *
 *  localData A pointer to the local data.  This pointer will be passed
 *     to the loader function.
 *
 *  loaderfcn A pointer to the function that will load data on the
 *     remote task into the send buffer.  This is where the caller
 *     gets to specify what fields need to be sent.
 *
 *  itemSize The number of bytes (per particle) that we want to transfer
 *
 *  remoteData An array that will be populated with the remote data
 *     for the particles in iLocal.  The caller is responsible to
 *     ensure there is memory allocated for itemSize*iLocalSize
 *     bytes.
 *
 *  The Loader Function
 *
 *  The caller must supply a loader function that copies the desired
 *  data fields into a comm buffer.  The arguments that will be passed
 *  to that function are:
 *
 *  commBuf The buffer where data should be placed.
 *
 *  nItems The number of particles for which data is to be loaded.
 *
 *  itemIndices The indices into the ddcSendList for particles for
 *      which data is requested
 *
 *  ddcSendList A segment of the ddc routing tables that lists the local
 *     particle indices that are sent to the task that is requesting data
 *
 *  In other words, if jj = ddcSendList[itemIndices[ii]], then jj is the
 *  local index (in STATE) of the iith particle for which data needs to
 *  be sent.  You can use jj to find the data that needs to be loaded
 *  into the buffer.
 *
 *  Implementation Notes:
 *
 *  This function works by taking advantage of the fact that remote
 *  particles are handled by the routing tables in ddc->CommInfo.  Hence
 *  for any remote particle that a given task knows about, it is
 *  possible to identify the task that owns that particle using
 *  ddc->CommInfo, send an request to that task, and then listen for the
 *  requested data to arrive.
 *
 *  We examine the list of request indices and determine which remote
 *  task the remote data is on, as well as how many requests will go to
 *  each remote task.  Since the local index will be meaningless to the
 *  remote task, we also find the position of the particle in the group
 *  of remote particles send from that task.
 *
 *  Next, we send the list of requests to each remote task.  There are
 *  three things to remember: 1.  We only need to communicate with tasks
 *  that sent us remote data so the list of tasks we communicate with is
 *  both limited and localized.  2.  Even if we are requesting no data
 *  from a task we still need to send the message because each task is
 *  expecting a message from each task it sends remote data to.  3.  The
 *  list of requests contains the index of the particle(s) within the
 *  block of remote particles sent from the remote task.
 *
 *  For each message a task receives, it uses MPI_Get_count to learn how
 *  many items are in the message.  If the answer is zero then there is
 *  nothing left to do since the message sender is not expecting a
 *  reply.  Otherwise, the loader function is called to load the
 *  apropriate data into a send buffer and the buffer is sent back to
 *  the task that originated the request.  The originating task receives
 *  the data into the appropriate slots of the array provided by the
 *  caller, and we're done.
 *
 *  The nSend and nRecv arrays should be at least as big as
 *  ddc->nRemoteProcessors.  They will be filled, respectively, with the
 *  number of particles for which this task sent requests for data to
 *  the ith remote processor or received requests for data from the ith
 *  remote processor.  This information is passed back to the caller
 *  since it can facilitate further communication about the remote data
 *  that is exchanged (such as sending back altered remote data).
 *
 *  Although this routine does strive to minimize the amount of
 *  communication, there is still non-trivial comm overhead.  Every time
 *  this routine is called each task sends a message to all tasks from
 *  which it recevies remote particles.  In the case that remote data is
 *  needed only rarely, these messages will consist overwhelmingly of
 *  messages with no data---a lot of traffic just to communicate that
 *  there is nothing to do.  
 *
 *  Usually we post recvs before sends to avoid unexpected msgs, but in
 *  this case the msgs are very small and there is likely to be enough
 *  load imbalance that there is no way to ensure all recvs are posted
 *  before the sends get started anyway.  The logic of the routine is
 *  easier if we do the sends first.
 *
 */
void getRemoteData(unsigned* iLocal, unsigned iLocalSize, void* localData,
 		   fcn_type loaderFcn,
		   unsigned itemSize, void* remoteData,
		   int* nSend, int* nRecv)
{
   assert(sorted(iLocal, iLocalSize));
   
   const int maxRequests = 512;
   static const int itemTag = 13;
   static const int dataTag = 14;

   unsigned itemSendBuf[iLocalSize];
   unsigned itemRecvBuf[maxRequests];
   char dataSendBuf[maxRequests*itemSize];
   
   DDC* ddc = getddc();

   MPI_Request itemSendReq[ddc->nRemoteProcessors];
   MPI_Request dataSendReq[ddc->nRemoteProcessors];


   unsigned kk=0;
   for (int ii=0; ii<ddc->nRemoteProcessors; ++ii)
   {
      unsigned recvStart = ddc->CommInfo[ii].RecvStart;
      unsigned recvEnd = recvStart + ddc->CommInfo[ii].nRecv;
      nSend[ii] = 0;
      for (; kk<iLocalSize && iLocal[kk]<recvEnd; ++kk)
      {
	 itemSendBuf[kk] = iLocal[kk] - recvStart;
	 ++nSend[ii];
      }
   }

   // post sends
   unsigned bufPos = 0;
   for (int ii=0; ii<ddc->nRemoteProcessors; ++ii)
   {
      unsigned* sendPtr = itemSendBuf + bufPos;
      unsigned dest = ddc->CommInfo[ii].domain_id;
      MPI_Isend(sendPtr, nSend[ii], MPI_UNSIGNED, dest, itemTag, COMM_LOCAL, itemSendReq+ii);
      bufPos += nSend[ii];
      assert(nSend[ii] < maxRequests);
   }

   unsigned nDataMsgs = 0; // number of remote tasks asking for non-zero amt of data.
   bufPos = 0;
   for (int ii=0; ii<ddc->nRemoteProcessors; ++ii)
   {
      MPI_Status status;
      unsigned id = ddc->CommInfo[ii].domain_id;
      MPI_Recv(itemRecvBuf, maxRequests, MPI_UNSIGNED, id, itemTag, COMM_LOCAL, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, nRecv+ii);
      if (nRecv[ii] == 0)
	 continue;
      int* ddcSendList = ddc->sendList + ddc->CommInfo[ii].SendStart;
      char* sendPtr = dataSendBuf + bufPos;
      unsigned nData = itemSize*nRecv[ii];
      loaderFcn(sendPtr, nRecv[ii], itemRecvBuf, ddcSendList, localData);
      MPI_Isend(sendPtr, nData, MPI_CHAR, id, dataTag, COMM_LOCAL, dataSendReq+nDataMsgs);
      ++nDataMsgs;
      bufPos += nData;
   }
   
   bufPos = 0;
   for (int ii=0; ii<ddc->nRemoteProcessors; ++ii)
   {
      if (nSend[ii] == 0)
	 continue;
      unsigned msgSize = nSend[ii] * itemSize;
      unsigned id = ddc->CommInfo[ii].domain_id;
      char* recvPtr = ((char*) remoteData) + bufPos;
      MPI_Recv(recvPtr, msgSize, MPI_CHAR, id, dataTag, COMM_LOCAL, MPI_STATUS_IGNORE);
      bufPos += msgSize;
   }

   // We need to wait for the sends before we exit the routine since the
   // send buffers are on the stack.  This also avoids a slow resource
   // leak that occurs if you never wait for (or free) the requests for
   // Isends and Irecvs.
   MPI_Waitall(ddc->nRemoteProcessors, itemSendReq, MPI_STATUSES_IGNORE);
   MPI_Waitall(nDataMsgs, dataSendReq, MPI_STATUSES_IGNORE);
}

