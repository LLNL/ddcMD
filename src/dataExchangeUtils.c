#include "dataExchangeUtils.h"

#include <assert.h>

#include "ddc.h"
#include "ddcMalloc.h"
#include "simulate.h"

static void findHaloDomains(double rcut, unsigned** haloList, unsigned* nHalo);
static unsigned selectHaloAtoms(STATE* state, double rCut, unsigned dest, unsigned* sendMap);

/** Builds a DATA_EXCHANGE_PARMS that will guarantee to communicate data
 *  for all remote particles that are within rCut of any local
 *  particle.
 *
 *  ALGORITHM:
 *
 *  Start by calling findHaloDomains to find all domains that are close
 *  enough to the local domain that there is some chance of
 *  communication.  This gives us a list of destination tasks and the
 *  number of tasks we will send to.
 *
 *  The sendOffset array is formed by looping over the destination
 *  tasks, and calling selectHaloAtoms for each such task.  This returns
 *  the number of atoms to send (so that the offset array can be
 *  populated), and fills in the corresponding block of the sendMap.
 */
DATA_EXCHANGE_PARMS mkRemoteAtomDep(double rCut)
{
   SIMULATE* simulate = simulate_getSimulate(NULL);
   unsigned nLocal = simulate->system->nlocal;
   STATE* state = simulate->system->collection->state;
   
   DATA_EXCHANGE_PARMS dp;
   
   findHaloDomains(rCut, &dp.destTask, &dp.nSend);
   dp.comm = COMM_LOCAL;
   dp.sendOffset = ddcMalloc((dp.nSend+1)*sizeof(unsigned));
   dp.sendMap = ddcMalloc(nLocal*dp.nSend*sizeof(unsigned));

   dp.sendOffset[0] = 0;
   for (unsigned ii=0; ii<dp.nSend; ++ii)
      dp.sendOffset[ii+1] = dp.sendOffset[ii] +
         selectHaloAtoms(state, rCut, dp.destTask[ii], dp.sendMap+dp.sendOffset[ii]);
   
   dep_negotiate(&dp);

   dp.recvMap = NULL;
   dp.initHook = dep_defaultInitHook;
   dp.freeHook = dep_defaultFreeHook;
   return dp;
}


/** Finds all of the remote domains that this task might interact with
 *  for an interaction radius of rcut.  Domains do not interact with
 *  themselves.
 *
 *  Ranks are given with respect to whatever comm describes the ranks
 *  for the ddc object.
 *
 *  The domain_overlap calculation uses the domain radii.  These radii
 *  are recomputed each time ddcSendRecvTables is called.  Hence, they
 *  may be a bit out of date.  However, we're just going to live with
 *  the problem since the error is likely to be small and would appear
 *  only at the tail of the correlation function.
 */
void findHaloDomains(double rcut, unsigned** haloList, unsigned* nHalo)
{
   DDC* ddc = getddc();
   int myId = ddc->domain_id;
   unsigned listSize = 0;
   unsigned listCapacity = 26;
   unsigned* list = ddcMalloc(listCapacity * sizeof(unsigned));
   for (int ii=0; ii<ddc->size; ++ii)
   {
      if (ii==myId)
         continue;
      if (domainset_overlap(&ddc->domains, ii, rcut))
      {
         if (listCapacity == listSize)
         {
            listCapacity *= 2;
            list = ddcRealloc(list, listCapacity * sizeof(unsigned));
         }
         list[listSize] = ii;
         ++listSize;
      }
   }
   *nHalo = listSize;
   *haloList = list;
}

/** Add atoms that might interact with the destination domain to the
 *  sendMap.  Returns the number of possibly interacting atoms.  Whether
 *  an atom possibly interacts is determined by domain_possibleRemote().
 */
unsigned selectHaloAtoms(STATE* state, double rCut, unsigned dest, unsigned* sendMap)
{
   double* rx = state->rx;
   double* ry = state->ry;
   double* rz = state->rz;
   unsigned nLocal = state->nlocal;
   const DOMAINX destDomain = getddc()->domains.domains[dest];
   
   THREE_VECTOR r;
   unsigned cnt = 0;
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      VSET(r, rx[ii], ry[ii], rz[ii]);
      if (domain_possibleRemote(destDomain, r, rCut))
         sendMap[cnt++] = ii;
   }
   return cnt;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
