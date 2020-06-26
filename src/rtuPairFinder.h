#ifndef RTU_AIR_FINDER_H
#define RTU_AIR_FINDER_H

#include "pairFinderGeom.h"
#include "simulate.h"

/** Combines a pairFinder with an array of atomic coordinates, rr to
 *  create a complete setup that is "ready to use" (hence rtu).  Whenever
 *  pfg is destroyed, rr must be removed from the heap. Use
 *  rtupf_destroy.  */

typedef struct rtu_pair_finder_st
{
   unsigned rBlk;
   unsigned nAtoms;
   double rcut;
   THREE_VECTOR* rr;
   PAIR_FINDER_GEOM pfg;
} RTU_PAIR_FINDER;

typedef struct displacement_struct
{
   double x,y,z,r2;
} DISPLACEMENT;

unsigned findNbrs(unsigned nPairs, RTU_PAIR_FINDER pairFinder,
		  unsigned iAtom, DISPLACEMENT* displacements);

RTU_PAIR_FINDER rtupf_create(SIMULATE* simulate, double rCut);

void rtupf_destroy(RTU_PAIR_FINDER* this);

int displacement_sort(DISPLACEMENT* d1, DISPLACEMENT* d2); 

#endif // ifndef RTU_AIR_FINDER_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
