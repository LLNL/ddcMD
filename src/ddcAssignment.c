#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <assert.h>
#include "ddc.h"
#include "ptiming.h"
#include "particle.h"
#include "heap.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include "mpiTypes.h"
#include "expandbuffer.h"
#include "domain.h"
#include "GridAssignmentObject.h"
//#include "back_communicate.h"
#include "alltoall_sparse.h"
#include "box.h"
#include "bioGid.h"

#include "simulate.h"
#include "bisectionCalc.h"

#include "preduce.h"
#include "units.h"
void voronoiCalcParticleDestinations(DDC* ddc, int verbose);

void ddcUpdateinfo(DDC*ddc);


static void ddcExchangeParticles(DDC* ddc, int verbose);
static void ddcSetOverLapID(DDC* ddc);

extern FILE *ddcfile;

/** sort_localDomain needs to be set to the local domain rank before
 *  calling sortParticlesByDomainId.  This allows particles in the local
 *  domain to be sorted to the front of the list regardless of the
 *  domain rank. */
static int sort_localDomain;
static int sortParticlesByDomainId(const PARTICLE* p1, const PARTICLE* p2);
int sortParticlesByGid(const PARTICLE* p1, const PARTICLE* p2);
static void ddcAffineUpdateCenters(DDC* ddc);


static int findCommInfo(int rank, DDC* ddc);


int ddcParticleDiffState(DDC*ddc,STATE *state) 
{
   int flag =0; 
   printf("%d: nlocal %d %d\n",getRank(0),state->nlocal,ddc->number_local);
   if ((int)state->nlocal != (int)ddc->number_local) flag =1; 
   if ((int)state->nion != (int)ddc->number_particles) flag+=2; 
   for (int i=0; i<state->nlocal;i++)
   {
      if (state->label[i] != ddc->particles[i].global_index)
      {
         flag+=4;
         break;
      }
   }
   return flag;
}
void ddcAssignment(DDC* ddc)
{
   static int first=1;
   ddc->Assignment_Called = 1; 

   if (first) timestampBarrier("Start ddcAssignment", COMM_LOCAL); 
   profile(ASSIGNMENT, START);

   //SYSTEM *sys = system_getSystem(NULL); 
   //STATE  *state = sys->collection->state; 
   ddc->number_particles = ddc->number_local;
   ddc->number_remote = 0;
   if (ddc->size ==0)
   {
      first =0;
      profile(ASSIGNMENT, END);
      if (first) timestampBarrier("End ddcAssignment", COMM_LOCAL); 
      return; 
   }
   ddcAffineUpdateCenters(ddc);
   profile(ASSIGNMENT_B0, START);
   if (first) timestampBarrier("Start domain_particle", COMM_LOCAL); 
   switch(ddc->loadBalance->itype)
   {
      case BISECTION:
         bisectionCalcParticleDestinations(ddc);
         break; 
      default:
         voronoiCalcParticleDestinations(ddc, first);
   }
   if (first) timestampBarrier("End domain_particle", COMM_LOCAL); 
   profile(ASSIGNMENT_B0, END);
   if (ddc->ddcRule != NULL) 
   {
      ddc->ddcRule->eval_ddcRule(ddc);
   }
   ddcExchangeParticles(ddc, first);
   profile(ASSIGNMENT_B0, END);
   if (first) timestampBarrier("End ddcAssignment", COMM_LOCAL); 
   profile(ASSIGNMENT, END);
   first = 0; 
}

/** This routine computes particle destinations by computing the minimum
 *  distance to any center in the list of domains that interact with
 *  this task (i.e., ddc->OverLapID).  Since this list is typically
 *  small, this "brute force" calculation is sufficiently fast.  If the
 *  OverLapIDs aren't set (probably because ddcSendRecvTables has never
 *  been called) then we need to bootstrap.  GridAssignmentObject will
 *  find destinations that are mostly correct (it ignores periodic
 *  boundary conditions).  Even if the destination is on the wrong side
 *  of a pbc, the particle will still be assigned to a domain that is
 *  close (actually adjacent I think) to the correct domain with pbc
 *  considered.  When the particles are sent to those domains the
 *  generated OverLapIDs will consist of a small list of domains.  This
 *  ensures that even for large numbers of particles on large numbers of
 *  tasks the initial assignment will be fast.
 */
void voronoiCalcParticleDestinations(DDC* ddc, int verbose)
{
   PARTICLE* ddcParticle = ddc->particles;
   if (ddc->OverLapID == NULL) 
   {
      GRID_ASSIGNMENT_OBJECT* gao =
         gao_init(ddc->size, &(ddc->domains.domains[0].center), sizeof(DOMAINX));
      for (unsigned ii=0; ii<ddc->number_particles; ++ii)
         ddcParticle[ii].domain_id = gao_nearestCenter(gao, ddcParticle[ii].r);
      gao_destroy(gao);
      for (unsigned ii=0; ii<ddc->pinnedSize; ++ii)
         ddcParticle[ddc->pinned[ii]].domain_id = ddc->domain_id;
      ddcExchangeParticles(ddc, verbose);
      ddcSetOverLapID(ddc);
   }

   //assert(ddc->OverLapID != NULL);
   int n = ddc->nOverLap+1; 
   int list[n]; 

   list[0]  = ddc->domain_id; 
   for (int i=1;i<n;i++) list[i] = ddc->OverLapID[i-1]; 
   for (unsigned i = 0; i < ddc->number_particles; i++) 
   {
      ddcParticle[i].domain_id = domainset_particle(&ddc->domains, ddcParticle[i].r,n,list);
   }
   for (unsigned ii=0; ii<ddc->pinnedSize; ++ii)
      ddcParticle[ddc->pinned[ii]].domain_id = ddc->domain_id;
}

/* Used when the taget array is a has table, for the back communication, 
 * to find the destination processor of count 
 */
#include "integer.h"
integer ivec_get_first(const void *ptr)
{
   return ((const int *) ptr)[0];
}

// Commenting out TRAGET_IS_HASH is a temporary fix to resolve a
// problem with the meteor impact simulations. When TARGET_IS_HASH
// is not defined, the old MPI_Alltoall() strategy is used for tasks
// to find out how many atoms to receive from other processors. When
// TARGET_IS_HASH is defined a new algorithm is used which uses considerably
// less memory and is faster on large numbers of tasks. The old strategy
// is perfectly fine, except when getting to >~1M tasks where it may
// cause memory problems.

//#define TARGET_IS_HASH

#ifdef TARGET_IS_HASH
/* { // Integer hash table to eliminate world-size target array */
typedef struct {
   int tabsize,nalloc,nused;
   int (*keydata)[2];
   int *next;
   int *hashtable;
} inthash;

static int isprime(int x) {
   int k = 3;
   if((x&1) == 0) return 0;
   while(k*k <= x) {
      if(x % k == 0) return 0;
      k += 2;
   }
   return 1;
}
static int nextprime(int x) {
   while(!isprime(x)) x++;
   return x;
}
static void inthash_new(inthash *vh,int nitems) {
   int i;
   vh->tabsize = nextprime(nitems);
   vh->nalloc = 0;
   vh->nused = 0;
   vh->keydata = NULL;
   vh->next = NULL;
   vh->hashtable = (int *) malloc(sizeof(int) * vh->tabsize);  
   for(i = 0; i<vh->tabsize; i++)
      vh->hashtable[i] = -1;
}
static void inthash_delete(inthash *vh) {
   free(vh->hashtable);
   free(vh->keydata);
   free(vh->next);
   vh->tabsize = vh->nalloc = vh->nused = 0;
   vh->hashtable = NULL;
   vh->keydata = NULL;
   vh->next = NULL;
}
static void * inthash_lookup_insert(inthash *vh,int key,int initializer) {
   int idx = ((size_t) key) % vh->tabsize;
   int p = vh->hashtable[idx];
   while(p >= 0 && vh->keydata[p][0] != key)
      p = vh->next[p];
   if(p < 0) {
      if(vh->nused >= vh->nalloc) {
         int n = vh->nalloc * 2;
         if(n < 100) n = 100;
         vh->next = (int *) realloc(vh->next,sizeof(*(vh->next)) * n);
         vh->keydata = (int (*)[2]) realloc(vh->keydata,sizeof(*(vh->keydata)) * n);
         vh->nalloc = n;
      }
      p = vh->nused++;
      vh->keydata[p][0] = key;
      vh->next[p] = vh->hashtable[idx];
      vh->hashtable[idx] = p;
      vh->keydata[p][1] = initializer;
   }
   return vh->keydata[p];
}
/*	} // End of integer hash table stuff */
#endif

/** This routine assumes that ddc->particles[ii].domain_id has been set
 *  to the rank of the task where particle ii is to be sent. */
void ddcExchangeParticles(DDC* ddc, int verbose)
{
   const int statTag=34;
   const int msgTag=35;

   profile(ASSIGNMENT_B2, START);
   if (verbose) timestampBarrier("Send Particle to Correct Domain", COMM_LOCAL); 
   PARTICLE* ddcParticle = ddc->particles;
   for (unsigned ii=0; ii<ddc->number_particles; ++ii)
      ddcParticle[ii].ifirst = ii;

   sort_localDomain = ddc->domain_id;
   qsort(ddcParticle, ddc->number_particles, sizeof(PARTICLE),                // Sort by destination 
         (int(*)(const void*,const void*))sortParticlesByDomainId);

#ifdef TARGET_IS_HASH
   inthash targethash;
   inthash_new(&targethash,200);
#else
   int target[ddc->size];
   for (int ii=0; ii<ddc->size; ++ii)
      target[ii] = 0;
#endif

   int startCapacity = 1024;
   int* start = ddcMalloc(startCapacity * sizeof(int));
   start[0] = 0;
   int lastid = ddc->domain_id;
   int nIsend = 0;
   for (unsigned i = 0; i < ddc->number_particles; i++)
   {
      if (ddcParticle[i].domain_id != lastid)
      {
         start[nIsend++] = i;
#ifdef TARGET_IS_HASH
         {
            int *target_ptr =
               inthash_lookup_insert(&targethash,
                     ddcParticle[i].domain_id,0 /* Initializer */);
            //assert(target_ptr[1] == target[p[i].domain_id]);
            target_ptr[1 /* 1 is data part, 0 is key part */] += 1;
         }
#else
         target[ddcParticle[i].domain_id] += 1;
#endif
         if (nIsend == startCapacity) // too small for next entry
         {
            startCapacity *= 2;
            start = ddcRealloc(start, startCapacity*sizeof(int));
         }
      }
      lastid = ddcParticle[i].domain_id;
   }
   assert(nIsend < startCapacity); // Ensure start is large enough. 
   start[nIsend] = ddc->number_particles; // start holds nIsend+1 entries

   // This check was put in during refactoring.  I can't imagine how it
   // could fail, but I put it in because I'm paranoid.
   if (ddc->number_particles == 0)
   {
      assert(nIsend == 0);
      assert(start[0] == 0);
   }

   profile(ASSIGNMENT_B2, END);
   profile(ASSIGNMENT_B3, START);
   profile(B3_1, START);

#ifdef TARGET_IS_HASH
   int nIrecv;
   {
      int maxcount;
      MPI_Allreduce(&targethash.nused,&maxcount,1,MPI_INT,MPI_MAX,COMM_LOCAL);
      maxcount *= 3;
      if(maxcount < 1000) maxcount = 1000;
      targethash.keydata = (int (*)[2]) realloc(targethash.keydata,
            sizeof(*(targethash.keydata)) * maxcount);
      targethash.nalloc = maxcount;

      {
         /* Communicate data to targets */
         if(0 && getRank(0) == 0) {
            printf("In %s() at %s:%d: -- Calling alltoall_sparse() to count receives...\n",
                  __func__,__FILE__,__LINE__);
            fflush(stdout);
         }
         int nrem = alltoall_sparse(0,ddc->size,ddc->domain_id,1,COMM_LOCAL,
               targethash.nused,
               targethash.nalloc,
               sizeof(*(targethash.keydata)),
               targethash.keydata,
               ivec_get_first);

         if(0 && getRank(0) == 0) {
            printf("In %s() at %s:%d: -- alltoall_sparse() complete.\n",
                  __func__,__FILE__,__LINE__);
            fflush(stdout);
         }


         /* Sum received items */
         nIrecv = 0;
         for(int i = 0; i<nrem; i++)
            nIrecv += targethash.keydata[i][1];
      }

      /* Don't need has any more */
      inthash_delete(&targethash);
   }
#else
#ifdef __bgq__
   MPI_Allreduce(MPI_IN_PLACE, target, ddc->size, MPI_INT, MPI_SUM, COMM_LOCAL);
#else
   /* If we don't have MPI_IN_PLACE, make temporary source array. */
   {
      int tmpsource[ddc->size];
      for(int k = 0; k<ddc->size; k++)
         tmpsource[k] = target[k];
      MPI_Allreduce(tmpsource, target, ddc->size, MPI_INT, MPI_SUM, COMM_LOCAL);
   }
#endif
   int nIrecv = target[getRank(0)];
   //assert(nIrecv == XnIrecv);
#endif
   profile(B3_1, END);

   static MPI_Request *recvRequest = NULL;
   static MPI_Request *sendRequest = NULL;
   // Each task will send nIsend messages and receive nIrecv messages.
   recvRequest = (MPI_Request*)ExpandBuffers(
         recvRequest, sizeof(MPI_Request), nIrecv, 512, LOCATION("ddcExchangeParticles"),"recvRequest");
   sendRequest = (MPI_Request*)ExpandBuffers(
         sendRequest, sizeof(MPI_Request), nIsend, 512, LOCATION("ddcExchangeParticles"),"sendRequest");

   // First communicate the message statistics: how many particles to receive and from where.
   if (verbose) timestampBarrier("Negotiating message routing", COMM_LOCAL);
   profile(B3_2, START);
   int* msgStat = ddcMalloc(2*nIrecv*sizeof(int));
   int* statBuf = ddcMalloc(2*nIsend*sizeof(int));
   for (int ii=0; ii<nIrecv; ++ii)
   {
      MPI_Irecv(msgStat+2*ii, 2, MPI_INT, MPI_ANY_SOURCE, statTag, COMM_LOCAL, recvRequest+ii);
   }
   for (int ii=0; ii<nIsend; ++ii)
   {
      int cnt = start[ii + 1] - start[ii];
      int dest = ddcParticle[start[ii]].domain_id;
      statBuf[2*ii] = cnt;
      statBuf[2*ii+1] = getRank(0);
      MPI_Isend(statBuf+2*ii, 2, MPI_INT, dest, statTag, COMM_LOCAL, sendRequest+ii);
   }

   int sendbufSize = start[nIsend] - start[0]; // total number of particles to send
   unsigned  sendbuf_blk;
   void* sendbuf = (void *)heapGet(&sendbuf_blk);
   heapEndBlock(sendbuf_blk, sendbufSize*sizeof(PARTICLE));
   PARTICLE* ParticlePtr = sendbuf;
   memcpy(sendbuf, ddcParticle + start[0], sendbufSize*sizeof(PARTICLE));

   ddc->nRemoteProcessors = nIsend; 
   ddc->CommInfo = (COMMINFO*) ExpandBuffers((void *)ddc->CommInfo, sizeof(COMMINFO), ddc->nRemoteProcessors, 128, LOCATION("ddcAssignment"),"ddc->CommInfo");
   for (int i=0; i<ddc->nRemoteProcessors; i++) 
   {
      ddc->CommInfo[i].domain_id=-1;
      ddc->CommInfo[i].nSend=0;
      ddc->CommInfo[i].SendStart=0;
      ddc->CommInfo[i].nRecv=0;
      ddc->CommInfo[i].RecvStart=0;
   }

   if (particleSizeinfo()>0) particleSortinfo((char *)&(ddcParticle->ifirst),sizeof(PARTICLE),ddc->number_local); 

   MPI_Waitall(nIsend, sendRequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nIrecv, recvRequest, MPI_STATUSES_IGNORE);
   profile(B3_2, END);

   if (verbose) timestampBarrier("Communicating particle data", COMM_LOCAL);
   profile(B3_3, START);

   ddcParticle=ddc->particles;
   int index = start[0];
   for (int ii=0; ii< nIrecv; ++ii)
   {
      int cnt = msgStat[2*ii];
      int source = msgStat[2*ii+1];
      MPI_Irecv(ddcParticle + index, cnt, particle_MPIType(), source, msgTag, COMM_LOCAL, recvRequest+ii);
      index += cnt;
   }

   index=0;
   for (int i = 0; i < nIsend; i++)
   {
      int cnt = start[i + 1] - start[i];
      int dest= ParticlePtr[index].domain_id;
      MPI_Isend(ParticlePtr + index, cnt, particle_MPIType(), dest, msgTag, COMM_LOCAL, sendRequest + i);
      ddc->CommInfo[i].domain_id=dest;
      ddc->CommInfo[i].nSend=cnt;
      ddc->CommInfo[i].SendStart=start[i];
      index += cnt;
   }

   index = start[0];
   for (int ii=0; ii<nIrecv; ++ii)
   {
      int cnt = msgStat[2*ii];
      int source = msgStat[2*ii+1];
      int iComm = findCommInfo(source, ddc);;
      ddc->CommInfo[iComm].nRecv=cnt;
      ddc->CommInfo[iComm].RecvStart=index;
      index += cnt;
   }
   ddcFree(start);
   profile(B3_3, END);

   profile(B3_4, START);
   MPI_Waitall(nIsend, sendRequest, MPI_STATUSES_IGNORE);
   MPI_Waitall(nIrecv, recvRequest, MPI_STATUSES_IGNORE);
   profile(B3_4, END);
   heapFree(sendbuf_blk); 

   ddc->number_particles = index;
   ddcFree(msgStat);
   ddcFree(statBuf);
   if (verbose) timestampBarrier("End  Receiving Particles", COMM_LOCAL); 
   profile(ASSIGNMENT_B3, END);
   profile(ASSIGNMENT_B5, START);
   ddc->number_local = ddc->number_particles;
   ddc->number_remote = 0;
   /* Sort particles by particle index to ensure that particle order is  independent of order that particles are received. */ 
   if (verbose) timestampBarrier("Start qsort", COMM_LOCAL); 
   for (unsigned i = 0; i < ddc->number_local; i++) ddcParticle[i].ifirst = i; /*Use ifirst as temp storage for the sort order*/
   qsort(ddcParticle, ddc->number_local, sizeof(PARTICLE), (int(*)(const void*,const void*))sortParticlesByGid); /*Sorts ddc copy of particle info by globalindex */
   if (verbose) timestampBarrier("End qsort", COMM_LOCAL); 
   if (verbose) timestampBarrier("Start ddcUpdateinfo", COMM_LOCAL); 
   if (particleSizeinfo()>0)
   {
      ddcUpdateinfo(ddc);
      particleSortinfo((char *)&(ddcParticle->ifirst),sizeof(PARTICLE),ddc->number_local);  /*Sorts MD particles structures in the same order as the ddc copy of particles i.e. sort using ifirst as the map*/
   }
   if (verbose) timestampBarrier("End ddcUpdateinfo", COMM_LOCAL); 
   profile(ASSIGNMENT_B5, END);
}

/** This code is similar to the first part of ddcSendRecvTables however,
 * it does not actually consider which particles interact with nbr
 * domains.  It also makes no attempt to avoid the MPI_Allgather that is
 * needed for all tasks to know the radii of all domains. */
void ddcSetOverLapID(DDC* ddc)
{
   int myId = ddc->domain_id;

   THREE_VECTOR center = domainset_getLocalCenter(&ddc->domains );
   ddc_getDomainSize(ddc, center);
   domainset_allGather(&ddc->domains);

   assert(ddc->rcut > 0);
   int nRemote = 0;
   double rcut_local = ddc->rcut;
   double dia = VSQ(box_get_diagonal(NULL));

   while (nRemote == 0) 
   {
      for (int ii=0; ii<ddc->size; ++ii)
      {
         if (ii == myId) continue;
         if ( domainset_overlap(&ddc->domains, ii, rcut_local) )
         {
            ddc->OverLapID = ddcRealloc(ddc->OverLapID, sizeof(int)*(nRemote+1));
            ddc->OverLapID[nRemote] = ii;
            nRemote++;
         }
      }
      if (rcut_local >  dia) break; 
      rcut_local = rcut_local * 1.26; /* 2**(1.0/3.0)) */
   }


   ddc->nOverLap = nRemote; 
}

/** Don't forget to set the static variable sort_localDomain before
 * calling this sort routine. */
int sortParticlesByDomainId(const PARTICLE* p1, const PARTICLE*p2)
{
   if (p1->domain_id == p2->domain_id) return 0;
   if (p1->domain_id == sort_localDomain) return -1;
   if (p2->domain_id == sort_localDomain) return 1;
   if (p1->domain_id > p2->domain_id) return 1;
   if (p1->domain_id < p2->domain_id) return -1;
   return 0;
}

int sortParticlesByGid(const PARTICLE*p1, const PARTICLE*p2)
{
   if (p1->global_index > p2->global_index) return 1;
   if (p1->global_index < p2->global_index) return -1;
   return 0;
}

/** This function returns the index in ddc->CommInfo where the recv
 *  information for the given rank can be stored.  It will either find
 *  the index in CommInfo where send information for the rank is already
 *  stored, or, if nothing is sent to the rank, it will expand the array
 *  and create a new entry in the array.
 *
 *  We go to this trouble because in ddcAssignment, the comm tables for
 *  send and recv are not symmetric.  For any rank to which particles
 *  will be sent, particles may or may not be received.  We may also
 *  receive particles from ranks to which none are sent.
 */
int findCommInfo(int rank, DDC* ddc)
{
   for (int ii=0; ii<ddc->nRemoteProcessors; ++ii)
      if (ddc->CommInfo[ii].domain_id == rank)
         return ii;
   int here = ddc->nRemoteProcessors++;
   ddc->CommInfo = (COMMINFO *)
      ExpandBuffers((void *)ddc->CommInfo,sizeof(COMMINFO),
            ddc->nRemoteProcessors, 128, LOCATION("ddcAssignment"),
            "ddc->CommInfo"); 
   ddc->CommInfo[here].domain_id=rank;
   ddc->CommInfo[here].nSend=0;
   ddc->CommInfo[here].SendStart=0;
   return here;	
}

void ddcAffineUpdateCenters(DDC* ddc)
{
   THREE_MATRIX hfac =
      box_getAffineUpdateHfac(NULL, ddc->centersAffineUpdateIndex);

   if ( matrix_equal_tol(hfac, I_3x3, 1e-14) )
      return;

   DOMAINX* d = ddc->domains.domains;
   for (int ii=0; ii<ddc->size; ++ii)
      d[ii].center = matrix_vector(hfac, d[ii].center);
}


/*
   --------------------------------
   A Brief History of ddcAssignment
   --------------------------------

   In the earliest versions of ddcMD, domain_size and domain_close were
   used to determine how many domains would be checked by domain_particle.
   For example, this method was in use in r74.  Presumeably a list of
   domains was also generated as a side effect of domain_close.  I haven't
   worked thorugh all of the forensics, but at some point this was
   abandoned in favor of attempts at lattice based methods, or using the
   existing ddc->OverlapID list.  Apparently the domain_size calculation
   was inadvertantly left in the code in spite of the fact that it was no
   longer used.  The domain_size calculation was finally removed in the
   major refactor in r1548.


   In the first version of ddcAssignment, once destinations were
   calculated, all particle messages were sent out using MPI_Isend.  The
   code then used an Iprobe strategy to discover which tasks wanted to
   exchange particle information.  This code worked by looking for messages
   from any task using MPI_Iprobe and receiving (with MPI_Recv) any that
   were present.  Every task would continue to probe and receive until no
   messages were found pending.  At this point, an all reduce was performed
   to see if all particles were accounted for.  If not, the code would
   usleep(1000) and go back into the Iprobe/recv loop.  This cycle would
   repeat up to ten times.

   There are three problems with this original method: 

   1.  First, the sends are posted before the receives.  This leads to a
   large number of unexpected messages which will cause some MPI
   implementation to allocate memory that will never be given back to the
   application.

   2.  There is no guarantee that the method will succeed.  There is no way
   to be certain that all of the messages will reach their targets in the
   10 possible probe loops.  Speaking practically, it always worked and it
   is hard to imagine how we could get past 10 usleeps and 10 global syncs
   (the allreduce) without delivering all of the Isends, but is is
   conceptually possible.

   3.  It is slow compared to a simpler implmentation.  In the replacement
   implementation we perform only one allreduce to determine how many
   messages are sent to each task.  Now we can recv exactly that many.
   This allows the recvs to post first and avoids the possibility multiple
   usleeps and allreduces.


   In r470 (15-Oct-2007) ddcAssignment was substantially rewritten.  This
   revision introduced a new method of particle exchange routing was based
   on an allreduce that would let each task know exactly how many tasks it
   would receive messages from.  Once this was known, tasks exchanged
   information to describe how many particles would come from each task.
   This allows each task to set up receive buffers for all particles that
   are sent and post recvs before sends.  It took a few revisions to thrash
   out all of the minor bugs, but this refactor proved very successful.

   In r490 (24-Oct-2007) the ddcAssignment function was moved from ddc.c to
   a new file ddcAssignment.c.  Memory allocation was changed in that the
   ddc->particles array and the send buffer sendbuf were no longer managed
   by ExpandBuffers but were taken from the scratch heap.  The old Iprobe
   based code was deleted.

   In r929 (5-Jun-2009) scratch heap safety improvements forced changes in
   the way heap allocations were handled.

   In r1548 (7-Apr-2011) ddcAssignment was again seriously refactored.
   This refactor split the destination calculation and data exchange tasks
   into two separate functions.  It also introduced a grid based method to
   bootstrap initial assignment.  This allows initial assignment to be fast
for large numbers of particles and tasks even on non-lattice domain
structures.  The previous brute force search over all centers for all
particles could take several minutes for large BGL or BGP runs.  This
was very painful when trying to debug a large scale problem.  The code
was also simplified by removing cruft that had accumulated over previous
revisions.

----------
Algorithm:
----------

As of r1548 this is the basic algorithm for ddcAssignment:

If there is only one task in the job, none of this is needed and we
effectively bail out immediately.  Otherwise,

            Part A:  Calculate destination domain for each particle:
            ddcCalculateParticleDestinations

            A1.  domain_set is called (unless a Hack flag is set).
            - The purpose of this call is to adjust the position of the 
            domain centers according to changes in box size or box shape. 
            - This is, as far as I know, the only mechanism currently
            implemented through which the box shape can change.

            A2.  The destination domain for each particle is now calculated and
            stored in the PARTICLE array ddc->particles.
            - The heavy lifting is performed by domainset_particle.  This 
            function calculates the minimum distance from the particle
            to a list of centers.  If ddc->OverLapID is non-null then we have
            a good list of candidate centers.  If not, (for example on the
                  initial assignment of atoms to tasks) we have to bootstrap the
            list.  
            - The bootstrapping procedure is to use GridAssignmentObject to 
            find domains for particles.  Particles are sent to domains found 
            by GridAssignmentObject by calling ddcExchangeParticles (the 
                  usual method of moving particles).  Because of simplifying
            assumptions made in GridAssignmentObject, particles may not 
            be assigned to the task in whose Voronoi cell they lie.  However,
            they will be (at worst) in a domain that is close.  This means
            we can use a bounding sphere approximation to build the
            ddc->OverLapID array (by calling ddcSetOverLapID) and be assured
            that the list of overlapping domains will be small.

            Part B:  Send particles to destination: ddcExchangeParticles

            B1.  Store an auxilliary sort key in ddcParticle[i].ifirst and sort the particles
            by destination domain.  The sort places local particles before
            remote particles.

            B2.  Loop over the particles.  Store the starting index for each block of
            particles headed for the same domain.  Also increment the target
            array for each task that will receive a message from this task.
            Perform an allreduce on the target array so that each task can know
            how many messages to expect.

            B3. Post recvs, sends, and waits to exchange message routing
            information.  Each task will learn which tasks will send data and
            how many particles.  
   - Before posting MPI_Waitall, copy the particle data into the send
buffer (take advantage of overlapping comm and compute)
   - Also call particleSortInfo which leverages the auxilliary sort
   key set previously to sort registered particle data that is 
   not included in the set of data directly stored in the
   PARTICLE structure.

   B4.  Post recvs, then sends for particle data.  As sends are posted,
   populate ddc->CommInfo.  Wait for particle data.

   B5.  Sort recived particles by gid.  Save an auxilliary sort key so that
   we can also sort any corresponding registered info.

   B6.  Call ddcUpdateinfo to exchange any registered particle
   information.  This uses the comm table created in ddc->CommInfo. 

   B7. Sort corresponding registered info.


   --------------------
   Performance counters
   --------------------

   ASSIGNMENT profiles the entire ddcAssignment routine
   ASSIGNMENT_B0 profiles ddcCalculateParticleDestinations
   DOMAIN_SET is reported as a sub-component of ASSIGNMENT_B0
ASSIGNMENT_B1 is now deleted (used to profile allreduce of nglobal)
   ASSIGNMENT_B2 covers sorting of particles by domainId and construction
   of target and start arrays
   ASSIGNMENT_B3 covers allreduce of target and all send/recv operations
   as well as associated Waitalls and computation of
   routing tables
   B3_1 covers allreduce of target
   B3_2 covers send/recv/wait of message size/rank info.  Also covers
   particleSortinfo.
   B3_3 covers sendrev of particle data
   B3_4 covers Waitall for particle data
   ASSIGNMENT_B4 was deleted in r470 (used to cover Waitall for 
         particle Isends)
   ASSIGNMENT_B5 covers sorting particles by gid, ddcUpdateinfo, and 
   particleSortinfo.



   ------
   Notes:
   ------

   1.  We should replace the call to domain_set with an alternate function
   that updates particle positions according to the change in the h-matrix
   since the last update.  This will allow us to handle changing box shape
   for any arrangement of centers.  (Presently only lattice arrangements
         work.)

   2.  In the r1548 refactor, the start, target, and target2 arrays were
   changed from compile time allocations of fixed size MAXPROC to run time
   allocations.  This probably moves these arrays from the system heap to
   the stack.  Since these arrays can be quite large it is worth asking if
   this was a good thing or if some alternate allocation method should be
   used.  

   3.  Is there any (portable) way to convert the allreduce of target to an
   in-place operation?  As I recall the in-place stuff doesn't work on blue
   gene.

   4.  Take a look at the (many) profile markers that are present in this
   code.  See which still make sense, which are reported or otherwise
   used. Reorganize as necessary.

   5.  The allreduce of the target array moves more data than is
   necessary.  A reduce scatter would be sufficient.  Is there a reasonable
   way to improve performance and memory use or is it really just a wash?

   6.  Why is there a 1000 microsecond sleep in ddcUpdateTables?
   (It has apparently been there, undisturbed, since r13.

    7.  It is worth asking whether ddcAssignment is really the right place
    to move the domain centers according to the box shape change.  Is there
    any chance that with a moving box simulation we could fail for a long
    time to trigger an update?

    8.  How is ddc->number_local initially set and how is it kept up to
    date?  The fact that ddc->number_remote is set to zero in ddcAssignment
    is a bit spooky.  For that matter, the fact that we have completely
    clobbered the comm tables is interesting.  It isn't entirely clear how
    the handshaking is done to let other ddc routines know that the comm
    tables are crap.

    9.  Why are the MPI_Request arrays static and managed by ExpandBuffers?
    We could easily put them on the stack, the scratch heap, or use
    malloc/free.  What advantage is the persistent assignment?  For first
    assignment in large simulations there are potentially a large number of
    requests. 

    10.  Think about whether ddcSetOverLapID should be more of a general
   purpose function that shares functionality with ddcSendRecvTables
   */


   /* Local Variables: */
   /* tab-width: 3 */
   /* End: */
