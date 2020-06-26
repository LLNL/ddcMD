#include "rtuPairFinder.h"
#include "system.h"
#include "heap.h"
#include "dataExchange.h"
#include "dataExchangeUtils.h"
#include "pairIterator.h"
#include "pairFinderGeom.h"
#include "preduce.h"

/** Find the atoms nearest to a central atom, iAtom, and populate a
 *  sorted list of their displacements relative to the central atom.
 *  Return the number of neighbors found.
 *
 *  This routine is required to find at least 4*nPairs neighbors.
 *  If there are not enough real neighbors, fill the remaining slots
 *  with synthetic data that has zero displacement, but an impossibly
 *  large distance.
 *
 *  It is assumed that the displacements array is large enough to hold
 *  as many atoms as the RTU_PAIR_FINDER will encounter.  Depending on
 *  rcut, this could be a large number of atoms.  The caller can
 *  check the returned number of neighbors to make sure there was no
 *  memory overrun.
 */
unsigned findNbrs(unsigned nPairs, RTU_PAIR_FINDER rtupf,
		  unsigned iAtom, DISPLACEMENT* displacements)
{
   unsigned nNbrs = 0;
   for (PAIR_ITERATOR iter=pfg_newIter(iAtom, &(rtupf.pfg));
	iter.ii!=iter.jj; pairIter_advance(&iter))
   {
      DISPLACEMENT* di = displacements+nNbrs++; 
      di->x = iter.rx;
      di->y = iter.ry;
      di->z = iter.rz;
      di->r2 = iter.r2;
   }
   
   qsort(displacements, nNbrs, sizeof(DISPLACEMENT),
	 (int (*)(const void *, const void *))displacement_sort);
   
   // if we didn't find enough real data, synthesize required amount.
   unsigned minNbrs = 4*nPairs;
   if (nNbrs < minNbrs)
   {
      for (unsigned ii=nNbrs; ii<minNbrs; ++ii)
      {
	 displacements[ii].x = 0.0; 
	 displacements[ii].y = 0.0; 
	 displacements[ii].z = 0.0; 
	 displacements[ii].r2 = 4 * rtupf.rcut * rtupf.rcut;
      }
      nNbrs = minNbrs;
   }
   return nNbrs;
}



/**
 *  This is analogous to a C++ constructor function.              
 *
 *  Combines a pairFinder with an array of atomic coordinates, rr to
 *  create a complete setup that is "ready to use" (hence rtu).
 *  Construct all needed data from scratch.                        
 */
RTU_PAIR_FINDER rtupf_create(SIMULATE* simulate, double rCut)
{
   SYSTEM* sys = simulate->system;
   STATE* state = sys->collection->state;
   unsigned nLocal = sys->nlocal;
   
   // Create the struct for data exchange (halo exchange). 
   DATA_EXCHANGE_PARMS dxParms = mkRemoteAtomDep(rCut);

   unsigned nRemote = dep_nRemote(dxParms);

   RTU_PAIR_FINDER rtupf;
   rtupf.nAtoms = nLocal+nRemote;

   // Create space for coordinates for all atoms in halo on the scratch
   // heap and copy all of the local coordinates into that space, then
   // execute the halo exchange of nremote coordinates.  Each task sends
   // local coordinates from ptr rr and receives remote coordinates to
   // rr+nlocal.
   rtupf.rr = heapGet(&(rtupf.rBlk));
   heapEndBlock(rtupf.rBlk, rtupf.nAtoms*sizeof(THREE_VECTOR));
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      rtupf.rr[ii].x = state->rx[ii];
      rtupf.rr[ii].y = state->ry[ii];
      rtupf.rr[ii].z = state->rz[ii];
   }
   dxParms.width = sizeof(THREE_VECTOR);
   dep_exchange(rtupf.rr, rtupf.rr+nLocal, &dxParms);

   dep_free(&dxParms);
   
   // We used to to make all the particles spatially compact.  This is
   // no longer necessary since the Geom takes care of it.

   rtupf.rcut = rCut;
   unsigned domainId = simulate->ddc->domain_id;
   const THREE_VECTOR* domainCenterP = &(simulate->ddc->domains.domains[domainId].center);
   PARTICLESET* ps = NULL;
   ps = ParticleSet(
      ps, &rtupf.rr[0].x, &rtupf.rr[0].y, &rtupf.rr[0].z, sizeof(THREE_VECTOR),
      NULL, 0, NULL, 0, &rtupf.nAtoms, &sys->nlocal, domainCenterP);

   rtupf.pfg = pfg_create(ps, nearestImage_fast, rtupf.rcut, box_getBox(NULL));

   return rtupf;
}


/** This is analogous to a C++ destructor function. 
 *  It must be used to prevent a memory leak.       
 */ 
void rtupf_destroy(RTU_PAIR_FINDER* this)
{
   pfg_destroy(&(this->pfg));
   heapFree(this->rBlk);
}


int displacement_sort( DISPLACEMENT *d1, DISPLACEMENT *d2)
{
   if  ( d1->r2 < d2->r2 ) return -1; 
   if  ( d1->r2 > d2->r2 ) return  1; 
   return  0; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
