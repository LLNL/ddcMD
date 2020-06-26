// should the create return type PAIR_FINDER_MODIFIED* instead?

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "object.h"
#include "pio.h"
#include "io.h"
#include "simulate.h"
#include "system.h"
#include "ddcMalloc.h"
#include "heap.h"
#include "crc32.h"
#include "units.h"
#include "preduce.h"

#include "dataExchange.h"
#include "dataExchangeUtils.h"
#include "pairIterator.h"
#include "pairFinder.h"
#include "pairFinderModified.h"

int getRank(int);

void pfm_destroy(PAIR_FINDER_MODIFIED * pfm)
{
 /////////////////////////////////////////////////////
 // This is analogous to a C++ destructor function. //
 // It must be used to prevent a memory leak.       //
 /////////////////////////////////////////////////////
   heapFree(pfm->iBlk);
   heapFree(pfm->vBlk);
   heapFree(pfm->rBlk);
   pfs_destroy(&(pfm->pfs));
   heapFree(pfm->pBlk);
}

PAIR_FINDER_MODIFIED * pfm_create(SYSTEM *sys, double rcut)
{
 ////////////////////////////////////////////////////////////////////
 // This is analogous to a C++ constructor function.               //
 ////////////////////////////////////////////////////////////////////
 // Group together Dave's simple pair finder struct with a list of //
 // the required atom positions, relative to an arbitrary origin.  //
 // Construct all needed data from scratch.                        //
 ////////////////////////////////////////////////////////////////////

// PAIR_FINDER_MODIFIED * pfm = ddcMalloc(sizeof(PAIR_FINDER_MODIFIED));
   unsigned int pBlkTmp;
   PAIR_FINDER_MODIFIED * pfm = heapGet(&pBlkTmp);
   heapEndBlock(pBlkTmp, sizeof(PAIR_FINDER_MODIFIED));
   pfm->pBlk=pBlkTmp;


 STATE * state = system_getState(NULL);
 unsigned nLocal = sys->nlocal;

 //////////////////////////////////////////////////////////
 // Create the struct for data exchange (halo exchange). //
 //////////////////////////////////////////////////////////

   DATA_EXCHANGE_PARMS dxParms = mkRemoteAtomDep(rcut);

   unsigned nRemote = dep_nRemote(dxParms);

   pfm->nAtoms = nLocal+nRemote;

// printf("ddt %d: nAtoms=%d\n", getRank(0), pfm->nAtoms);


 ////////////////////////////////////////////////////
 // Create space for nion coordinates in the heap. //
 ////////////////////////////////////////////////////
 pfm->rr = heapGet(&(pfm->rBlk));
 heapEndBlock(pfm->rBlk, pfm->nAtoms*sizeof(THREE_VECTOR));

 ///////////////////////////////////////////////////
 // Create space for nion velocities in the heap. //
 ///////////////////////////////////////////////////
 pfm->vr = heapGet(&(pfm->vBlk));
 heapEndBlock(pfm->vBlk, pfm->nAtoms*sizeof(THREE_VECTOR));

  
 ////////////////////////////////////////////////
 // Create space for nion species in the heap. //
 ////////////////////////////////////////////////
 pfm->atomType = heapGet(&(pfm->iBlk));
 heapEndBlock(pfm->iBlk, pfm->nAtoms*sizeof(int));

 /////////////////////////////////////////////////////////////
 // Copy all of the local coordinates into that heap space. //
 /////////////////////////////////////////////////////////////
 for (unsigned ii=0; ii<nLocal; ++ii)
 {
    pfm->rr[ii].x = state->rx[ii];
    pfm->rr[ii].y = state->ry[ii];
    pfm->rr[ii].z = state->rz[ii];
    pfm->vr[ii].x = state->vx[ii];
    pfm->vr[ii].y = state->vy[ii];
    pfm->vr[ii].z = state->vz[ii];
    pfm->atomType[ii] = state->atomtype[ii];
 }

 ////////////////////////////////////////////////////////
 // Execute the halo exchange of nremote coordinates   //
 // Each task sends local coordinates from ptr rr      //
 // and receives remote coordinates to rr+nlocal.      //
 ////////////////////////////////////////////////////////
 dxParms.width = sizeof(THREE_VECTOR);
 dep_exchange(pfm->rr, pfm->rr+nLocal, &dxParms);
 dep_exchange(pfm->vr, pfm->vr+nLocal, &dxParms);

 /////////////////////////////////////////////
 // Repeat halo exchange for species types. //
 /////////////////////////////////////////////
 dxParms.width = sizeof(int);
 dep_exchange(pfm->atomType, pfm->atomType+nLocal, &dxParms);

 /////////////////////////////////////////////////////
 // make atoms spatially compact;                   //
 // Define a local origin, r0, and map              //
 // relative positions to lie inside one unit cell. //
 /////////////////////////////////////////////////////
 THREE_VECTOR r0 = pfm->rr[0];
 for (unsigned ii=0; ii<pfm->nAtoms; ++ii)
 {
    pfm->rr[ii].x -= r0.x;
    pfm->rr[ii].y -= r0.y;
    pfm->rr[ii].z -= r0.z;
    nearestImage_fast(&(pfm->rr[ii].x), &(pfm->rr[ii].y), &(pfm->rr[ii].z));
 }

 /////////////////////////////////////////
 // Create struct for simple pair list. //
 /////////////////////////////////////////
 pfm->rcut  = rcut;
 pfm->rcut2 = rcut*rcut;
 pfm->pfs = pfs_create(pfm->rr, pfm->nAtoms, nearestImage_fast, pfm->rcut);

   dep_free(&dxParms);

 //////////////////////////////////////
 // Return the simple pair list      //
 // along with all atom coordinates. //
 //////////////////////////////////////
 return pfm;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
