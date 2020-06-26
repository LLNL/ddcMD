#include "ddc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "object.h"
#include "particle.h"
#include "box.h"
#include "heap.h"
#include "ddcMalloc.h"
#include "expandbuffer.h"
#include "md2ddc.h"

#include "system.h"

#include "rectimer.h"
#include "functimer.h"

#define TAG 1


static MPI_Request* _remoterequest = NULL;
static MPI_Request* _Sremoterequest = NULL;


void (*forceaccum[])(double *force, int nlist, int *list) ={forceaccum0,forceaccum1,forceaccum2,forceaccum3,forceaccum4,forceaccum5,forceaccum6,forceaccum7};

void expandStaticBuffers(const int comm_info_size)
{
   _remoterequest  = (MPI_Request *)ExpandBuffers(_remoterequest,  sizeof(MPI_Request), comm_info_size, 512, LOCATION("ddcUpdate"), "remoterequest");
   _Sremoterequest = (MPI_Request *)ExpandBuffers(_Sremoterequest, sizeof(MPI_Request), comm_info_size, 512, LOCATION("ddcUpdate"), "Sremoterequest");
}

void ddcUpdate(DDC*ddc)
{
   COMMINFO* comm_info      = ddc->CommInfo;
   const int comm_info_size = ddc->nRemoteProcessors;

   int nRecv;
   void *RecvPtr;
   int SendStart, nlist, remoteID, RecvStart;
   MPI_Datatype RecvType;
   int *list;
   profile(UPDATE, START);
   unsigned send_blk; 
   double* send = (double *)heapGet(&send_blk);
   unsigned sendSize = 0;
   expandStaticBuffers(comm_info_size);
	int nR=0; 
   for (int i = 0; i < comm_info_size; i++)
   {
	  
      nRecv = comm_info[i].nRecv;
	  if (nRecv==0) continue; 
      remoteID = comm_info[i].domain_id;
      RecvStart = comm_info[i].RecvStart;
      comm_info[i].RecvPtr=RecvPtr=GetRecvLocation(nRecv,RecvStart,&(comm_info[i].RecvType));
      RecvType = comm_info[i].RecvType ;
      MPI_Irecv(RecvPtr, 1, RecvType, remoteID, TAG, COMM_LOCAL, _remoterequest + nR++ );
   }
   int n=3; 
   int nS=0; 
   for (int i = 0; i < comm_info_size; i++)
   {
      nlist = comm_info[i].nSend;
	  if (nlist == 0) continue ;
      remoteID = comm_info[i].domain_id;
      SendStart = comm_info[i].SendStart;
      list = ddc->sendList + SendStart;
      fillsendbuf(send, nlist, list);
      MPI_Isend(send, n*nlist, MPI_DOUBLE, remoteID, TAG, COMM_LOCAL, _Sremoterequest + nS++);
      send += n*nlist;
      sendSize += n*nlist;
   }
   heapEndBlock(send_blk, sendSize*sizeof(double));
   
   profile(UPDATEWAIT, START);
   MPI_Waitall(nS, _Sremoterequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nR, _remoterequest, MPI_STATUSES_IGNORE);
   profile(UPDATEWAIT, END);
   heapFree(send_blk); 
   profile(UPDATE, END);
}
void ddcUpdateVelocity(DDC*ddc)
{
   COMMINFO* comm_info=ddc->CommInfo;
   const int comm_info_size = ddc->nRemoteProcessors;

   int nRecv;
   void *RecvPtr;
   int SendStart, nlist, remoteID, RecvStart;
   MPI_Datatype RecvType;
   int *list;
   profile(UPDATE, START);
   unsigned send_blk; 
   double* send = (double *)heapGet(&send_blk);
   unsigned sendSize=0;
   expandStaticBuffers(comm_info_size);
	int nR=0; 
   for (int i = 0; i < comm_info_size; i++)
   {
      nRecv = comm_info[i].nRecv;
	  if (nRecv ==0) continue; 
      remoteID = comm_info[i].domain_id;
      RecvStart = comm_info[i].RecvStart;
        comm_info[i].RecvPtr=RecvPtr=GetVelocityLocation(nRecv,RecvStart,&(comm_info[i].RecvType));
      RecvType = comm_info[i].RecvType ;
      MPI_Irecv(RecvPtr, 1, RecvType, remoteID, TAG, COMM_LOCAL, _remoterequest + nR++);
   }
   int n=3; 
	int nS=0; 
   for (int i = 0; i < comm_info_size; i++)
   {
      nlist = comm_info[i].nSend;
	  if (nlist ==0) continue; 
      remoteID = comm_info[i].domain_id;
      SendStart = comm_info[i].SendStart;
      list = ddc->sendList + SendStart;
      Velocityfillsendbuf(send, nlist, list);
      MPI_Isend(send, n*nlist, MPI_DOUBLE, remoteID, TAG, COMM_LOCAL, _Sremoterequest + nS++);
      send += n*nlist;
      sendSize += n*nlist;
   }
   heapEndBlock(send_blk, sendSize*sizeof(double));

   profile(UPDATEWAIT, START);
   MPI_Waitall(nS, _Sremoterequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nR, _remoterequest,  MPI_STATUSES_IGNORE);
   profile(UPDATEWAIT, END);
   heapFree(send_blk); 
   profile(UPDATE, END);
}

void ddcUpdateForce(DDC*ddc, unsigned calculate)
{
   if (ddc->PairMode == 0) return;

	STARTTIMER;
	static thandle
	  t_postrecv = NULL,
	  t_postsend = NULL,
	  t_waitrecv = NULL,
	  t_waitsend = NULL,
	  t_forceacc = NULL;

   COMMINFO* comm_info=ddc->CommInfo;
   const int comm_info_size = ddc->nRemoteProcessors;

   int n, size, cnt;
   int nlist, remoteID, Start;
   int *list;
   profile(P_UPDATEFORCE, START);
   cnt =0; 
   for (int i = 0; i < comm_info_size; i++) cnt += comm_info[i].nSend;
   int ndouble=3;
   if ((calculate & pPotentialEnergyMask) > 0)  ndouble += 1; 
   if ((calculate & pVirialMask) > 0)  ndouble += 1; 
   if ((calculate & pSionMask) > 0)  ndouble += 6; 
   
   //int ndouble=11;
   size = ndouble*sizeof(double); 
   unsigned recvbuf_blk; 
   void* recvbuf=heapGet(&recvbuf_blk);
   heapEndBlock(recvbuf_blk, cnt*size);
   double* bufPtr = (double *)recvbuf;
   expandStaticBuffers(comm_info_size);
	int nR=0; 

	t_postrecv = rectimer_start(t_postrecv,"f.updt.recv.post");
   for (int i = 0; i < comm_info_size; i++)
   {
      nlist = comm_info[i].nSend;
	  if (nlist ==0) continue; 
      remoteID = comm_info[i].domain_id;
      MPI_Irecv(bufPtr, ndouble*nlist, MPI_DOUBLE, remoteID, TAG, COMM_LOCAL, _remoterequest+nR++); 
      bufPtr += ndouble*nlist;
   }
	rectimer_stop(t_postrecv);

	t_postsend = rectimer_start(t_postsend,"f.updt.send.post");
	int nS=0; 
   for (int i = 0; i < comm_info_size; i++)
   {
      n = comm_info[i].nRecv;
	  if (n == 0) continue; 
      remoteID = comm_info[i].domain_id;
      Start = comm_info[i].RecvStart;
      void* Ptr = GetForceLocation(n, Start, &(comm_info[i].ForceType),calculate);
      MPI_Isend(Ptr, 1, comm_info[i].ForceType, remoteID, TAG, COMM_LOCAL, _Sremoterequest + nS++);
   }
	rectimer_stop(t_postsend);

	t_waitsend = rectimer_start(t_waitsend,"f.updt.send.wait");
   MPI_Waitall(nS, _Sremoterequest, MPI_STATUSES_IGNORE);
	rectimer_stop(t_waitsend);

	t_waitrecv = rectimer_start(t_waitrecv,"f.updt.recv.wait");
   MPI_Waitall(nR, _remoterequest, MPI_STATUSES_IGNORE);
	rectimer_stop(t_waitrecv);

	t_forceacc = rectimer_start(t_forceacc,"f.updt.accum");
   bufPtr = (double *)recvbuf;
   for (int i = 0; i < comm_info_size; i++)
   {
      nlist = comm_info[i].nSend;
      Start = comm_info[i].SendStart;
      list = ddc->sendList + Start;
      forceaccum[calculate](bufPtr, nlist, list);
      bufPtr += ndouble*nlist;
   }
	rectimer_stop(t_forceacc);

   heapFree(recvbuf_blk); 
   profile(P_UPDATEFORCE, END);

	STOPTIMER;
}
void ddcUpdateinfo(DDC*ddc)
{
   COMMINFO* comm_info=ddc->CommInfo;
   const int comm_info_size = ddc->nRemoteProcessors;

   int cnt, nRecv, size;
   void *RecvPtr;
   int SendStart, nlist, remoteID, RecvStart;
   MPI_Datatype RecvType;
   cnt =0; 
   for (int i = 0; i < comm_info_size; i++) cnt += comm_info[i].nSend;
   size = particleSizeinfo(); 
   unsigned sendbuf_blk; 
   char* sendbuf = NULL;
   char* send = sendbuf = heapGet(&sendbuf_blk); 
   expandStaticBuffers(comm_info_size);

   int nRecvRequest,nSendRequest; 
   nRecvRequest=nSendRequest=0; 
   send = sendbuf; 
   for (int i = 0; i < comm_info_size; i++)
   {
      nlist = comm_info[i].nSend;
      SendStart = comm_info[i].SendStart;
      particleFillBuffer(send,SendStart,nlist);  
      send += nlist*size;
   }
   heapEndBlock(sendbuf_blk, send-sendbuf);

	{
	  /*
		 Need to call resize() rather than particleAllocateinfo().
		 Otherwise the auxiliary particle arrays may be shrunk, while
		 on next resize() they will be assumed to be at least as large
		 as the lasrgest size previously seen by resize().
	  */
	  resize(ddc->number_local,2,system_getSystem(NULL)->collection->state);
	}

   for (int i = 0; i < comm_info_size; i++)
   {
      remoteID = comm_info[i].domain_id;
      nRecv = comm_info[i].nRecv;
      RecvStart = comm_info[i].RecvStart;
      if (nRecv > 0 ) 
      {
         comm_info[i].RecvPtr=RecvPtr=GetLocationinfo(nRecv,RecvStart,&(comm_info[i].RecvType));
         RecvType=comm_info[i].RecvType; 
         MPI_Irecv(RecvPtr,1,RecvType,remoteID,TAG,COMM_LOCAL,_remoterequest+nRecvRequest);   
         nRecvRequest++; 
      }
   }
   send = sendbuf; 
   for (int i = 0; i < comm_info_size; i++)
   {
      remoteID = comm_info[i].domain_id;
      nlist = comm_info[i].nSend;
      SendStart = comm_info[i].SendStart;
      if (nlist > 0 ) 
      {
         MPI_Isend(send, nlist*size, MPI_BYTE, remoteID, TAG, COMM_LOCAL, _Sremoterequest + nSendRequest); 
         nSendRequest++; 
         send += nlist*size;
      }
   }
   MPI_Waitall(nSendRequest, _Sremoterequest, MPI_STATUSES_IGNORE); 
   MPI_Waitall(nRecvRequest, _remoterequest, MPI_STATUSES_IGNORE); 
   heapFree(sendbuf_blk);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
