#ifndef PAIR_FINDER_H
#define PAIR_FINDER_H

#include "pairIterator.h"
#include "three_algebra.h"


typedef void (*nearestImageFcn) (double* x, double* y, double* z);

/** The collection of atoms given to the pfs_create function
 * should be compact. */

typedef struct PairFinderSimple_st
{
   unsigned nAtoms;
   THREE_VECTOR* r;
   THREE_VECTOR rMin;
   THREE_VECTOR rMax;
   THREE_VECTOR drInv;
   THREE_INT    nCell;
   
   double rCut2;
   int* firstAtomInCell;
   int* nextAtomInCell;
   nearestImageFcn nearestImage;
} PAIR_FINDER_SIMPLE;


PAIR_FINDER_SIMPLE pfs_create(THREE_VECTOR* r,
			      unsigned nAtoms,
			      nearestImageFcn nearestImage,
			      double rCut);
void pfs_destroy(PAIR_FINDER_SIMPLE* grid);
PAIR_ITERATOR pfs_newIter(int iAtom, PAIR_FINDER_SIMPLE* this);
void pfs_advanceIter(PAIR_ITERATOR* iter);


#endif // #ifndef PAIR_FINDER_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
