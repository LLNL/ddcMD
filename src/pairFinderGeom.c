#include "pairFinderGeom.h"
#include "ddcMalloc.h"

typedef struct PairFinderGeomState_st
{
   int stridedouble;
   double* xptr;
   double* yptr;
   double* zptr;
   GEOM_PINFO* pinfo;
   GEOMBOX* box0;
   int* nList;
   int nn;
   int iBox;
   double rCut2;
   nearestImageFcn nearestImage;
   THREE_VECTOR ri;
} PAIR_FINDER_GEOM_STATE;

/** The newly created PAIR_FINDER_GEOM takes over ownership of the
 *  PARTICLESET* that is passed to the constructor.  It will be
 *  destroyed by pfg_destroy.  
 */
PAIR_FINDER_GEOM pfg_create(PARTICLESET* ps, nearestImageFcn nearestImage, double rCut, BOX_STRUCT *box)
{
   PAIR_FINDER_GEOM pfg;
   pfg.geom = geom_init(rCut, box);
   pfg.ps = ps;
   pfg.nearestImage = nearestImage;
   GeomBox(pfg.geom, ps);
   return pfg;
}

void pfg_destroy(PAIR_FINDER_GEOM* this)
{
   geom_free(this->geom);
   ddcFree(this->ps);
}

PAIR_ITERATOR pfg_newIter(int iAtom, PAIR_FINDER_GEOM* this)
{
   PAIR_ITERATOR iter;
   iter.ii = iAtom;
   iter.advance = pfg_advanceIter;
   PAIR_FINDER_GEOM_STATE* state = ddcMalloc(sizeof(PAIR_FINDER_GEOM_STATE));
   iter.pairFinderState = state;

   state->stridedouble = this->ps->stridedouble;
   state->xptr  = this->ps->xptr;
   state->yptr  = this->ps->yptr;
   state->zptr  = this->ps->zptr;
   state->pinfo = this->geom->pinfo;
   state->box0  = this->geom->pinfo[iAtom].box;
   state->nList = state->box0->nlist;
   state->nn    = state->box0->nn;
   state->iBox  = 0;
   state->rCut2 = this->geom->rcut*this->geom->rcut;
   state->nearestImage = this->nearestImage;
   int stride = state->stridedouble;
   state->ri.x  = *(double*)((char*)(state->xptr) + iAtom*stride);
   state->ri.y  = *(double*)((char*)(state->yptr) + iAtom*stride);
   state->ri.z  = *(double*)((char*)(state->zptr) + iAtom*stride);

   iter.jj = (state->box0+state->nList[state->iBox])->first;
   while (iter.jj < 0)
   {
      ++state->iBox;
      if (state->iBox == state->nn)
      {
	 iter.jj = iter.ii;
	 ddcFree(state);
	 return iter;
      }
      GEOMBOX* box = state->box0 + state->nList[state->iBox];
      iter.jj = box->first;
   }

   double rjx  = *(double*)((char*)(state->xptr) + iter.jj*stride);
   double rjy  = *(double*)((char*)(state->yptr) + iter.jj*stride);
   double rjz  = *(double*)((char*)(state->zptr) + iter.jj*stride);

   iter.rx = state->ri.x - rjx;
   iter.ry = state->ri.y - rjy;
   iter.rz = state->ri.z - rjz;
   state->nearestImage(&iter.rx, &iter.ry, &iter.rz);
   iter.r2 = iter.rx*iter.rx + iter.ry*iter.ry + iter.rz*iter.rz;

   if (iter.ii == iter.jj || iter.r2 > state->rCut2)
      pfg_advanceIter(&iter);

   return iter;
}


void pfg_advanceIter(PAIR_ITERATOR* iter)
{
   PAIR_FINDER_GEOM_STATE* state = iter->pairFinderState;

   do
   {
      iter->jj = state->pinfo[iter->jj].next;
      while (iter->jj < 0)
      {
	 ++state->iBox;
	 if (state->iBox == state->nn)
	 {
	    iter->jj = iter->ii;
	    ddcFree(state);
	    return;
	 }
	 GEOMBOX* box = state->box0 + state->nList[state->iBox];
	 iter->jj = box->first;
      }
      int stride = state->stridedouble;
      double rjx  = *(double*)((char*)(state->xptr) + iter->jj*stride);
      double rjy  = *(double*)((char*)(state->yptr) + iter->jj*stride);
      double rjz  = *(double*)((char*)(state->zptr) + iter->jj*stride);
      
      iter->rx = state->ri.x - rjx;
      iter->ry = state->ri.y - rjy;
      iter->rz = state->ri.z - rjz;
      state->nearestImage(&iter->rx, &iter->ry, &iter->rz);
      iter->r2 = iter->rx*iter->rx + iter->ry*iter->ry + iter->rz*iter->rz;
      
   } while (iter->r2 > state->rCut2 || iter->ii == iter->jj);
}


