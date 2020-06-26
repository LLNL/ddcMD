#include "gidShuffle.h"
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include "random.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))



typedef struct GidShuffle_st
{
   double   order;
   gid_type label;
} GIDSHUFFLE;

static int gidShuffleLess(const void* a, const void* b)
{
   GIDSHUFFLE* aa = (GIDSHUFFLE*) a;
   GIDSHUFFLE* bb = (GIDSHUFFLE*) b;
   if (aa->order< bb->order) return -1;
   if (aa->order==bb->order) return 0;
   return 1;
}

static void fisherYatesShuffle64(gid_type* a, unsigned n)
{
   if (n == 0) return;
   for (unsigned ii=0; ii<n-1; ++ii)
   {
      unsigned jj = ii + (drand48() * (n - ii));
      gid_type tmp = a[ii];
      a[ii] = a[jj];
      a[jj] = tmp;
   }
}

static void fisherYatesShuffle32(unsigned* a, unsigned n)
{
   if (n == 0) return;
   for (unsigned ii=0; ii<n-1; ++ii)
   {
      unsigned jj = ii + (drand48() * (n - ii));
      unsigned tmp = a[ii];
      a[ii] = a[jj];
      a[jj] = tmp;
   }
}

/** Sorts unsigned data */
static int unsignedSortFunction(const void* av, const void* bv)
{
	const unsigned* a = (const unsigned*)av;
	const unsigned* b = (const unsigned*)bv;
	if (*a > *b) return 1;
	if (*a < *b) return -1;
	return 0;
}

static unsigned* buildDestArray(unsigned nItems, unsigned nSplits)
{
   unsigned myRank = getRank(0);
   unsigned nTasks = getSize(0);
   
   unsigned myTargets[nSplits];
   unsigned* tasks = ddcMalloc(nTasks*sizeof(unsigned));
   for (unsigned ii=0; ii<nTasks; ++ii)
      tasks[ii] = ii;
   
   for (unsigned ii=0; ii<nSplits; ++ii)
   {
      unsigned short seed[3];
      if (myRank == 0)
      {
	 seed[0] = drand48() * USHRT_MAX;
	 seed[1] = drand48() * USHRT_MAX;
	 seed[2] = drand48() * USHRT_MAX;
      }
      MPI_Bcast(seed, 3, MPI_UNSIGNED_SHORT, 0, COMM_LOCAL);
      seed48(seed);
      fisherYatesShuffle32(tasks, nTasks);
      myTargets[ii] = tasks[myRank];
   }
   ddcFree(tasks);

   qsort(myTargets, nSplits, sizeof(unsigned), unsignedSortFunction);

   unsigned* dest = ddcMalloc(nItems*sizeof(unsigned));
   unsigned blockSize = ceil(nItems / (nSplits*1.0));
   for (unsigned ii=0; ii<nSplits; ++ii)
      for (unsigned jj=0; jj<blockSize; ++jj)
	 if (ii*blockSize + jj < nItems)
	    dest[ii*blockSize+jj] = myTargets[ii];

   return dest;
}


/** This routine should shuffle the gids the same way regardless of how
 *  many tasks are used.  See collection_buildLattice in collection.c for
 *  more information.
 *
 *  This routine is typically extremely slow on more than a few thousand
 *  tasks.  For large numbers of tasks use gidShuffleLarge instead.
 *
 *  Algorithm:
 *
 *  #  For each local particle, associate the particle label with a
 *     random number in (0,1).  Call this number the order.
 *     - prand48 is used to make the result independent of nTasks.
 *  #  Globally sort the labels according to the number just picked.
 *     - The global sort works by associating each task with a segment
 *       in the range (0,1).  Labels with order in the segment are
 *       sent to the corresponding task using the assignArray
 *       utility.  This is what generates (potentially) large
 *       amounts of point to point traffic.  Now the labels
 *       are sorted by order locally.  The result is a virtual
 *       global array of labels that have been randomly
 *       reordered by the initial random number.
 *  #  Call distributeArray so that the local number of elements
 *       of the virtual global array is the same as the number
 *       of particles.  Then just reassign the labels from the
 *       virtual global array. 
 */
void gidShuffleSmall(COLLECTION* collection,  LONG64 seed)
{ 
   unsigned nLocal = collection->state->nlocal;
   gid_type* label = collection->state->label;

   timestampBarrier("Performing gid shuffle", COMM_LOCAL);

   unsigned nTasks = getSize(0);
   unsigned* count = ddcMalloc(nTasks*sizeof(unsigned));
   unsigned* countSum = ddcMalloc(nTasks*sizeof(unsigned));
   unsigned capacity = nLocal;
   unsigned* dest = ddcMalloc(nLocal*sizeof(unsigned));
   GIDSHUFFLE* newLabel = ddcMalloc(capacity*sizeof(GIDSHUFFLE));
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      PRAND48_STATE handle = prand48_init(label[ii], seed, 0xe1d2c3b4a5f7llu);
      for (unsigned jj=0; jj<1000; ++jj)
	 prand48(&handle);
      newLabel[ii].label = label[ii];
      newLabel[ii].order = prand48(&handle);
   }
   qsort(newLabel, nLocal, sizeof(GIDSHUFFLE), gidShuffleLess);

   for (unsigned ii=0; ii<nTasks; ++ii)
      count[ii] = 0;
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      dest[ii] = newLabel[ii].order*nTasks;
      count[dest[ii]] += 1;
   }
   
   MPI_Allreduce(count, countSum, nTasks, MPI_INT, MPI_SUM, COMM_LOCAL);

   capacity = MAX(nLocal, countSum[getRank(0)]);
   newLabel = ddcRealloc(newLabel, capacity*sizeof(GIDSHUFFLE));

   ddcFree(countSum);
   ddcFree(count);
   
   unsigned nHave = nLocal;
   assignArray((unsigned char*) newLabel, &nHave, capacity, sizeof(GIDSHUFFLE), dest, 1, COMM_LOCAL);
   ddcFree(dest);
   qsort(newLabel, nHave, sizeof(GIDSHUFFLE), gidShuffleLess);

   distributeArray((unsigned char*) newLabel, nHave, nLocal, sizeof(GIDSHUFFLE), COMM_LOCAL);

   for (unsigned ii=0; ii<nLocal; ++ii)
      label[ii] = newLabel[ii].label;
   
   ddcFree(newLabel);
   
   timestampBarrier("gid shuffle complete", COMM_LOCAL);

}

void gidShuffleLarge(COLLECTION* collection, LONG64 seed)
{ 
   const unsigned nSplits = 5;
   const unsigned nShuffles = 10;
   timestampBarrier("Performing gid shuffle", COMM_LOCAL);
   unsigned short seedVec[3];
   seedVec[0] = seed & 0xFFFF; seed = seed >> 16;
   seedVec[1] = seed & 0xFFFF; seed = seed >> 16;
   seedVec[2] = seed & 0xFFFF;
   seed48(seedVec);
   unsigned nLocal = collection->state->nlocal;
   gid_type* label = collection->state->label;

   unsigned maxCapacity = ceil(nLocal/(nSplits*1.0))*nSplits;
   unsigned send = maxCapacity;
   MPI_Allreduce(&send, &maxCapacity, 1, MPI_INT, MPI_MAX, COMM_LOCAL);
   gid_type* newLabel = ddcMalloc(maxCapacity*sizeof(gid_type));
   for (unsigned ii=0; ii<nLocal; ++ii)
      newLabel[ii] = label[ii];
   
   unsigned nHave = nLocal;
   for (unsigned iShuffle=0; iShuffle<nShuffles; ++iShuffle)
   {
      fisherYatesShuffle64(newLabel, nHave);
      unsigned* dest = buildDestArray(nHave, nSplits);
      assignArray((unsigned char*) newLabel, &nHave, maxCapacity, sizeof(gid_type), dest, 1, COMM_LOCAL);
      ddcFree(dest);
   }
   distributeArray((unsigned char*) newLabel, nHave, nLocal, sizeof(gid_type), COMM_LOCAL);
   
   for (unsigned ii=0; ii<nLocal; ++ii)
      label[ii] = newLabel[ii];
   
   ddcFree(newLabel);
   timestampBarrier("gid shuffle complete", COMM_LOCAL);
}

void gidShuffle(COLLECTION* collection, LONG64 seed)
{
   if (getSize(0) > 4096)
      gidShuffleLarge(collection, seed);
   else
      gidShuffleSmall(collection, seed);
}


void gidReset(COLLECTION* collection)
{
   unsigned nLocal = collection->state->nlocal;

   timestampBarrier("Performing gid reset", COMM_LOCAL);

   unsigned nTasks = getSize(0);
   unsigned myRank = getRank(0);
   unsigned* nAtoms = ddcMalloc(nTasks*sizeof(unsigned));

   MPI_Allgather(&nLocal, 1, MPI_INT, nAtoms, 1, MPI_INT, COMM_LOCAL);

   gid_type offset = 0;
   for (unsigned ii=0; ii<myRank; ++ii)
      offset += nAtoms[ii];

   for (unsigned ii=0; ii<nLocal; ++ii)
      collection->state->label[ii] = offset+ii;

   timestampBarrier("Finished gid reset", COMM_LOCAL);

   ddcFree(nAtoms);
}

