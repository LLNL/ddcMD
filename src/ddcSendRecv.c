#include "ddc.h"
#include <assert.h>
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
#include "preduce.h"
#include "mpiTypes.h"
#include "units.h"
#include "loadBalance.h"


#include "rectimer.h"
#include "functimer.h"


#define TAG 1


static int *_remoteSendStart = NULL;
static int *_remoteSendSize = NULL;
static int *_remoteRecvSize = NULL;

void ddcFindPairs(DDC*ddc);
int domain_possibleRemote_by_plane(DDC*ddc, THREE_VECTOR *u, THREE_VECTOR v, double *dist); 
double domain_separating_plane(THREE_VECTOR *c0, THREE_VECTOR *c1,THREE_VECTOR *u);
void *GetRecvLocation(int nRecv, int RecvStart, MPI_Datatype*RECVTYPE);

void ddcSendRecvTables(DDC*ddc)
{
  //if (getSize(0)==1){return;}
  
  STARTTIMER;
  static thandle t_comm = NULL,t_tx_1 = NULL,t_tx_2 = NULL,t_loop = NULL,t_dom,t_allg = NULL;
	int i, tag = 1;
	void *sendBuf = NULL;
	unsigned sendBuf_blk; 
   const int id = ddc->domain_id;

	profile(COMMLIST, START);
   profile(COMMDOMS, START);
   t_dom = rectimer_start(t_dom,"COMMDOMS");

   // get domain size from local particles positions
   { // limit scope
      THREE_VECTOR center = domainset_getLocalCenter( &ddc->domains );
      ddc_getDomainSize(ddc, center);
      //const double radius = ddc_getDomainSize(ddc, center);
      //domainset_setLocalRadius( &ddc->domains, radius );
   }

   int use_by_plane=0; 
   double minspan = box_get_minspan(NULL); 
   if (ddc->loadBalance->itype  != BISECTION) 
   {
      t_allg = rectimer_start(t_allg,"C-D : allgather");
      double max_radius; 
      if (ddc->itype==ARBITRARY ) 
      {
         domainset_allGather( &ddc->domains );
         max_radius = domainset_getMaxRadius( &ddc->domains );
      }
      else
      {
         max_radius = domainset_setRadiusToMax( &ddc->domains );
      }
      
      rectimer_stop(t_allg);

      // test that simulation box is big enough that domains can't see themselves.
      if (2.0*max_radius < 0.5*minspan*0.999 && ddc->Assignment_Called)
         use_by_plane=1; 
   }

   rectimer_stop(t_dom);
   profile(COMMDOMS, END);

	profile(COMMCALC, START);
    
	/* Can not use "by plane" method with bisection load balancer */
   //if (ddc->loadBalance->itype  == BISECTION) use_by_plane = 0;

   
   //  Find "interacting" domains
   DOMAINX *local_domain = ddc->domains.local_domain;
   if (ddc->loadBalance->itype  != BISECTION)
   {
      int remoteCapacity=0; 
      int nremote = 0;
      for (int idremote = 0 ; idremote < ddc->size; idremote++)
      {
         DOMAINX *remote_domain = ddc->domains.domains+idremote;
         if (idremote == id) continue;
         if (! domainset_overlap( &ddc->domains, idremote, ddc->rcut) ) continue;
         if ( ddc->comm_mode == TWOSIDED || domain_side(local_domain, remote_domain) > 0 ) 
         {
            if (nremote >= remoteCapacity) 
            {
               ddc->OverLapID = ddcRealloc(ddc->OverLapID, sizeof(int)*(nremote + 64));
               remoteCapacity=nremote + 64;
            }
            ddc->OverLapID[nremote] = idremote;
            nremote++;
         }
      }
      ddc->nOverLap = nremote; 
   }


   profile(COMMDOMS, START);
   int nremote = ddc->nOverLap;
   DOMAINX nbrDomains[nremote];
   {
      MPI_Request request[2*nremote];
      MPI_Datatype domainType = domainx_MPIType();
      for (int ii=0; ii<nremote; ++ii)
      {
         int remoteId = ddc->OverLapID[ii];
         MPI_Irecv(nbrDomains+ii, 1, domainType, remoteId, 0, COMM_LOCAL, request+ii);
         MPI_Isend(local_domain,   1, domainType, remoteId, 0, COMM_LOCAL, request+nremote+ii);
      }
      MPI_Waitall(2*nremote, request, MPI_STATUSES_IGNORE);
   }
   profile(COMMDOMS, END);



   _remoteSendStart = (int *)ExpandBuffers(_remoteSendStart, sizeof(int), nremote , 512, LOCATION("ddcSendRecvTables"),"_remoteSendStart");
   _remoteSendSize = (int *)ExpandBuffers(_remoteSendSize, sizeof(int), nremote , 512, LOCATION("ddcSendRecvTables"),"_remoteSendSize");
   _remoteRecvSize = (int *)ExpandBuffers(_remoteRecvSize, sizeof(int), nremote , 512,LOCATION("ddcSendRecvTables"), "_remoteRecvSize");

   profile(FIND_REMOTE, START);
   sendBuf = heapGet(&sendBuf_blk);
   unsigned sendBufSize = 0;
   unsigned sendBufCapacity = heapBlockSize(sendBuf_blk) / sizeof(PARTICLE);
   int cnt=0;
   int kk=0; 
   for (int kRemote=0;kRemote<nremote;kRemote++) 
   {
      THREE_VECTOR u[2]; 
      if (use_by_plane || ddc->forceUsePlanes)
      {
         DOMAINX *remote_domain = nbrDomains+kRemote; 
         double dist = domainset_separatingPlane(local_domain, remote_domain, u);
         if (dist >= 0.5*minspan*0.999) use_by_plane=0;
      }
      cnt = 0;
      for (unsigned j = 0; j < ddc->number_local; j++)
      {
         int flag; 
         double dummy; 

         if (use_by_plane || ddc->forceUsePlanes)
            flag=domain_possibleRemote_by_plane(ddc, u, ddc->particles[j].r,&dummy);
         else
            flag = domain_possibleRemote(nbrDomains[kRemote], ddc->particles[j].r, ddc->rcut);
         if (flag)
         {
            ++sendBufSize;
            assert(sendBufSize < sendBufCapacity);
            ((PARTICLE *) sendBuf)[kk + cnt] = ddc->particles[j];
            ((PARTICLE *) sendBuf)[kk + cnt].local_index = j; //dfr
            cnt++;
         }
      }
      _remoteRecvSize[kRemote] = -1;
      _remoteSendStart[kRemote] = kk;
      _remoteSendSize[kRemote] = cnt;
      kk += cnt;
   }
   heapEndBlock(sendBuf_blk, sendBufSize*sizeof(PARTICLE));
   profile(FIND_REMOTE, END);
   profile(COMMCALC, END);
   profile(COMMCOMM, START);
   MPI_Request remoterequest[nremote+1],Sremoterequest[nremote+1]; 
   t_comm = rectimer_start(t_comm,"COMMCOMM");
   t_tx_1 = rectimer_start(t_tx_1,"C-C : tx phase 1");

   for (i = 0; i < nremote; i++)
      MPI_Irecv(_remoteRecvSize + i, 1, MPI_INT, ddc->OverLapID[i], TAG, COMM_LOCAL, remoterequest + i);
   for (i = 0; i < nremote; i++)
      MPI_Isend(_remoteSendSize + i, 1, MPI_INT, ddc->OverLapID[i], TAG, COMM_LOCAL, Sremoterequest +i);
   MPI_Waitall(nremote, Sremoterequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nremote, remoterequest, MPI_STATUSES_IGNORE);

   rectimer_stop(t_tx_1);

   int number_remote =0; 
   for (i = 0; i < nremote; i++)
      number_remote += _remoteRecvSize[i];
   ddc->number_particles = ddc->number_local+number_remote;
   ddc->number_remote = number_remote;
   int k = ddc->number_local;

   t_tx_2 = rectimer_start(t_tx_2,"C-C : tx phase 2");
   for (i = 0; i < nremote; i++)
   {
      int remoteID = ddc->OverLapID[i];
      cnt = _remoteRecvSize[i];
      if (cnt > 0)
         MPI_Irecv(ddc->particles + k, cnt, particle_MPIType(), remoteID, tag, COMM_LOCAL, remoterequest + i);
      k += cnt;
   }
   for (i = 0; i < nremote; i++)
   {
      PARTICLE *ParticlePtr = ((PARTICLE *) sendBuf) + _remoteSendStart[i];
      cnt = _remoteSendSize[i];
      int remoteID = ddc->OverLapID[i];
      if (cnt > 0)
         MPI_Isend(ParticlePtr, cnt, particle_MPIType(), remoteID, tag, COMM_LOCAL, Sremoterequest + i);
   }
   MPI_Waitall(nremote, Sremoterequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nremote, remoterequest, MPI_STATUSES_IGNORE);
   rectimer_stop(t_tx_2);

   t_loop = rectimer_start(t_loop,"C-C : loop");

   int nRemoteAtoms =0; 
   if (nremote>0)
      nRemoteAtoms = _remoteSendSize[nremote-1] + _remoteSendStart[nremote-1]  - _remoteSendStart[0];
   ddc->sendList = (int *)ExpandBuffers((void *)ddc->sendList, sizeof(int), nRemoteAtoms, 4096, LOCATION("ddcSendRecvTables"),"ddc->sendList");

   ddc->nRemoteProcessors = nremote;
   ddc->nSend = nRemoteAtoms;
   // build send list and offsets.
   for (int ii=0; ii<nRemoteAtoms; ++ii)
      ddc->sendList[ii] = ((PARTICLE*)sendBuf)[_remoteSendStart[0]+ii].local_index;
   k = ddc->number_local;

   // setup ddc->CommInfo
   ddc->CommInfo = (COMMINFO *) ExpandBuffers((void *)ddc->CommInfo, sizeof(COMMINFO), nremote, 128,
         LOCATION("ddcSendRecvTables"),"ddc->CommInfo");
   for (int ii=0; ii<nremote; ++ii)
   {
      ddc->CommInfo[ii].domain_id = ddc->OverLapID[ii];
      ddc->CommInfo[ii].nSend     = _remoteSendSize[ii];
      ddc->CommInfo[ii].SendStart = _remoteSendStart[ii] - _remoteSendStart[0];
      cnt = _remoteRecvSize[ii];
      ddc->CommInfo[ii].nRecv = cnt;
      ddc->CommInfo[ii].RecvStart = k;
      ddc->CommInfo[ii].RecvPtr = GetRecvLocation(cnt, k, &(ddc->CommInfo[ii].RecvType));
      k += cnt;
   }
   
   THREE_VECTOR center = domainset_getLocalCenter(&ddc->domains);
   for (unsigned i=0;i<ddc->number_particles;i++)
   {
      THREE_VECTOR delta; 
      VOP2(delta,=,ddc->particles[i].r,-,center); 
      nearestImage_fast(&delta.x,&delta.y,&delta.z);
      VOP2(ddc->particles[i].r,=,delta,+,center); 
   }
   
   rectimer_stop(t_loop);
   rectimer_stop(t_comm);

   profile(COMMCOMM, END);
   profile(FINDPAIR, START);
   if (ddc->pruneSends == 1) ddcFindPairs(ddc);  
   profile(FINDPAIR, END);
   heapFree(sendBuf_blk); 
   profile(COMMLIST, END);

   STOPTIMER;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
