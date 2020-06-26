#include "ddcenergy.h"

#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>
#include "three_algebra.h"
#include "object.h"
#include "ddc.h"
#include "heap.h"
#include "particle.h"
#include "parityHandler.h"
#include "utilities.h"
#include "neighbor.h"
#include "ptiming.h"
#include "units.h"
#include "mpiUtils.h"

int sortParticlesByGid(const PARTICLE* p1, const PARTICLE* p2);
void ddcUpdate(DDC *ddc);
void constructList(SYSTEM *sys, double rcut); 

static PARTICLESET *_particleset = NULL;

PARTICLESET *getParticleSet()
{
      return _particleset; 
}
void setParticleSet(DDC* ddc, SYSTEM* sys)
{
   STATE* state = sys->collection->state;
   _particleset = ParticleSet(_particleset, state->rx, state->ry, state->rz, 
                              8, state->label, 8, state->atomtype, 4, &sys->nion,
                              &(ddc->number_local),
                              &(ddc->domains.domains[ddc->domain_id].center));
}
void zeroRpairs(SYSTEM*sys)
{
   for (unsigned ij = 0; ij < sys->neighbor->nRpair; ij++)
   {
      sys->neighbor->Rpairs[ij].fp = vzero;
      sys->neighbor->Rpairs[ij].e=0.0;
   }
}

int check4updateNeighbor(SYSTEM* sys)
{
   profile(VERLET, ACCUM);
   int update = neighborCheck(sys->neighbor, _particleset);
   profile(VERLET, END);

   //profile(MPIMISC, ACCUM);
   int updatesum;
   MPI_Allreduce(&update, &updatesum, 1, MPI_INT, MPI_SUM, COMM_LOCAL);
   //profile(MPIMISC, END);

   if (updatesum > 0) update = 1;

   return update;
}

int evalUpdateFlag(DDC* ddc, SYSTEM*sys, int ForcedUpdate)
{
   int update = 0;
   if (ForcedUpdate == -1) return update; 
   if (ddc->lastUpdate == NO_DDC_UPDATE || TEST0(sys->loop, ddc->updateRate) || ForcedUpdate) return 1;
	if (ddc->updateRate==0) update = check4updateNeighbor(sys);
   return update;
}

/** We have to be somewhat clever with our heap usage here.  We can't
 * know in advance how many particles ddcAssignment is going to send us
 * or how many remote particles we will get from ddcSendRecvTables.
 * However we also can't leave our heap block open since those two
 * subroutines also use the heap.  Hence the best we can do is make a
 * guess regarding how much space we will need and check to make sure we
 * didn't exceeed it.  */

void ddcUpdateTables(DDC*ddc, SYSTEM*system)
{
   profile(UPDATETABLES, START);
   if (getSize(0)  == 1) 
   {
      ddc->update = 0;
      if (ddc->lastUpdate ==NO_DDC_UPDATE) ddc->lastUpdate = system->loop;
      return; 
   }

	MPI_Bcast(&system->box->volume, 1, MPI_DOUBLE, 0, COMM_LOCAL);
	unsigned pblk;
	ddc->particles=heapGet(&pblk);
	size_t pblkSize = heapBlockSize(pblk)/2;
	heapEndBlock(pblk, pblkSize); 
   size_t needSize = (ddc->number_local) * sizeof(PARTICLE);
	if (needSize > pblkSize ) heapTooSmall(__FILE__, __LINE__, needSize ,pblkSize);
	ddcGetParticles(ddc);
   ddcAssignment(ddc);
   needSize = (ddc->number_local) * sizeof(PARTICLE);
	if (needSize > pblkSize ) heapTooSmall(__FILE__, __LINE__, needSize ,pblkSize);
 //return ; //  hack 
	system->neighbor->center = domainset_getLocalCenter( &ddc->domains );
	if (getSize(0) > 1) ddcSendRecvTables(ddc);
   needSize = (ddc->number_particles) * sizeof(PARTICLE);
	if (needSize > pblkSize ) heapTooSmall(__FILE__, __LINE__, needSize ,pblkSize);
	resize(ddc->number_local + ddc->number_remote, 2, system->collection->state);
	ddcPutParticles(ddc);
	heapFree(pblk); 
	usleep(1000);
	ddc->lastUpdate = system->loop;
	ddc->update = 0;
	system->nlocal = ddc->number_local;
	system->nremote = ddc->number_remote;
	system->nion = ddc->number_local + ddc->number_remote;
	if (ddc->domain_id == 0) printf("ddcUpdateTables @ loop=%"PRId64"\n", system->loop);

   profile(UPDATETABLES, END);
}
int  ddcUpdateAll(DDC*ddc, SYSTEM*sys, ETYPE *e, int ForcedUpdate)
{
   STATE *state = sys->collection->state;
   state->nlocal = sys->nlocal;
   state->nion = sys->nion;
   setParticleSet(ddc, sys);
   int update = evalUpdateFlag(ddc, sys, ForcedUpdate);
   if (update)
   {
      ddcUpdateTables(ddc, sys);
      profile(VERLET, START);
      state = sys->collection->state;
      state->nlocal = sys->nlocal;
      state->nion = sys->nion;
      setParticleSet(ddc, sys);
      sys->neighbor->lastUpdate=sys->loop; 
      if ((sys->neighborTableType & NEIGHBORTABLE_GPU)>0) 
      {
         constructList(sys, ddc->rcut); 
      }
      if ((sys->neighborTableType & NEIGHBORTABLE_FAT)> 0)
      { 
         neighborRef(sys->neighbor, _particleset);
         neighbors(sys->neighbor, _particleset); 
         zeroRpairs(sys); 
      }
      if ((sys->neighborTableType & NEIGHBORTABLE_SKINNY)>0)  
      {
         neighborRef(sys->neighbor, _particleset);
         neighbors1(sys->neighbor, _particleset); 
      }
      if ((sys->moleculeClass != NULL) && (ddc->lastUpdate == sys->loop)) moleculeScanState(sys->moleculeClass, state);
      profile(VERLET, END);
      return update; 
   }
if (getSize(0) > 0) ddcUpdate(ddc);  //For update steps the remote particle has been exchanged in the ddcUpdateTables but for the other steps ddcUpdate needs to be called.
   if ((sys->neighborTableType & NEIGHBORTABLE_FAT)>0)  
   {
      profile(UPDATEPAIR, START);
      pairUpdate(sys->neighbor, _particleset); 
      zeroRpairs(sys); 
      profile(UPDATEPAIR, END);
   }
   return update; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
