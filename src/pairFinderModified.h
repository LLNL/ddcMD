// should the create return type PAIR_FINDER_MODIFIED* instead?
#ifndef PAIR_FINDER_MODIFIED_H
#define PAIR_FINDER_MODIFIED_H

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

#include "dataExchange.h"
#include "dataExchangeUtils.h"
#include "pairIterator.h"
#include "pairFinder.h"

typedef struct pair_finder_modified_st
{
 ///////////////////////////////////////////////
 // Group Dave's simple pair finder together  //
 // with its array of atomic coordinates, rr, //
 // and species atomType.                     //
 // When pfs is destroyed, rr and atomType    //
 // must be removed from the heap.            //
 // Use pfm_destroy.                          //
 //////////////////////////////////////////////////////////
 // e.g., use as follows to copy x-separation of pairs   //
 // PAIR_FINDER_MODIFIED pfm;                            //
 // jx = 0;                                              //
 // for (PAIR_ITERATOR iter=pfs_newIter(ii, &(pfm.pfs)); //
 //     iter.ii!=iter.jj; pairIter_advance(&iter))       //
 //     x[jx++] = pfm.rr[iter.jj].x - pfm.rr[iter.ii].x; //
 //////////////////////////////////////////////////////////

 double rcut, rcut2;
 PAIR_FINDER_SIMPLE pfs;
 unsigned pBlk, rBlk, iBlk, vBlk, nAtoms;
 THREE_VECTOR * rr, * vr;
 int * atomType;
} PAIR_FINDER_MODIFIED;

void pfm_destroy(PAIR_FINDER_MODIFIED * pfm);

PAIR_FINDER_MODIFIED * pfm_create(SYSTEM *sys, double rcut);

#endif// PAIR_FINDER_MODIFIED_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
