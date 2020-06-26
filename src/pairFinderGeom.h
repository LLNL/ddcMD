#ifndef PAIR_FINDER_GEOM_H
#define PAIR_FINDER_GEOM_H

#include "pairIterator.h"
#include "box.h"
#include "geom.h"
#include "three_algebra.h"
#include "pairFinder.h" // nearestImageFcn typedef

/** This is an implementation of the PairFinder interface that uses the
 *  GEOM class to do all of the heavy lifting.  */
typedef struct PairFinderGeom_st
{
   GEOM* geom;
   PARTICLESET* ps;
   nearestImageFcn nearestImage;
} PAIR_FINDER_GEOM;

PAIR_FINDER_GEOM pfg_create(PARTICLESET* ps, nearestImageFcn nearestImage, double rCut, BOX_STRUCT *box);

void pfg_destroy(PAIR_FINDER_GEOM* this);
PAIR_ITERATOR pfg_newIter(int iAtom, PAIR_FINDER_GEOM* this);
void pfg_advanceIter(PAIR_ITERATOR* iter);


#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
