#include "pairFinder.h"

#include <assert.h>
#include "ddcMalloc.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))


typedef struct PairFinderSimpleState_st
{
   int iCell; // index to the 27 cells we search
   int nCells;
   int cellList[27];
   PAIR_FINDER_SIMPLE* grid;
} PAIR_FINDER_SIMPLE_STATE;

static int
pfs_mkCellList(unsigned iAtom, int* cellList, PAIR_FINDER_SIMPLE* this);


PAIR_FINDER_SIMPLE pfs_create(THREE_VECTOR* r,
			      unsigned nAtoms,
			      nearestImageFcn nearestImage,
			      double rCut)
{
   PAIR_FINDER_SIMPLE grid;
   grid.nAtoms = nAtoms;
   grid.r = r;
   grid.nearestImage = nearestImage;
   grid.rCut2 = rCut*rCut;
   
   grid.rMin = grid.rMax = r[0];
   for (unsigned ii=1; ii<nAtoms; ++ii)
   {
      grid.rMin.x = MIN(grid.rMin.x, r[ii].x);
      grid.rMax.x = MAX(grid.rMax.x, r[ii].x);
      grid.rMin.y = MIN(grid.rMin.y, r[ii].y);
      grid.rMax.y = MAX(grid.rMax.y, r[ii].y);
      grid.rMin.z = MIN(grid.rMin.z, r[ii].z);
      grid.rMax.z = MAX(grid.rMax.z, r[ii].z);
   }
   
   THREE_VECTOR boxSize;
   boxSize.x = grid.rMax.x - grid.rMin.x;
   boxSize.y = grid.rMax.y - grid.rMin.y;
   boxSize.z = grid.rMax.z - grid.rMin.z;
   
   assert(boxSize.x > 2*rCut &&
	  boxSize.y > 2*rCut &&
	  boxSize.z > 2*rCut );  // this should be an error msg.

   grid.nCell.x = boxSize.x / rCut;
   grid.nCell.y = boxSize.y / rCut;
   grid.nCell.z = boxSize.z / rCut;
   grid.drInv.x = grid.nCell.x / boxSize.x;
   grid.drInv.y = grid.nCell.y / boxSize.y;
   grid.drInv.z = grid.nCell.z / boxSize.z;

   unsigned numberCells = grid.nCell.x * grid.nCell.y * grid.nCell.z;
   grid.firstAtomInCell= ddcMalloc(sizeof(double)*numberCells);
   grid.nextAtomInCell = ddcMalloc(sizeof(int)*grid.nAtoms);
   for (unsigned ii=0; ii<numberCells; ++ii)
      grid.firstAtomInCell[ii] = -1;
   for (unsigned ii=0; ii<nAtoms; ++ii)
      grid.nextAtomInCell[ii] = -2;
   
   for (unsigned ii=0; ii<nAtoms; ++ii)
   {
      int ix = (r[ii].x - grid.rMin.x) * grid.drInv.x;
      int iy = (r[ii].y - grid.rMin.y) * grid.drInv.y;
      int iz = (r[ii].z - grid.rMin.z) * grid.drInv.z;
      while (ix>=grid.nCell.x) ix -= grid.nCell.x;
      while (iy>=grid.nCell.y) iy -= grid.nCell.y;
      while (iz>=grid.nCell.z) iz -= grid.nCell.z;
      unsigned index = ix + grid.nCell.x*(iy + grid.nCell.y*iz);
      assert(index<numberCells);
      grid.nextAtomInCell[ii] = grid.firstAtomInCell[index];
      grid.firstAtomInCell[index] = ii;
   }
   return grid;
}


void pfs_destroy(PAIR_FINDER_SIMPLE* grid)
{
   ddcFree(grid->nextAtomInCell); 
   ddcFree(grid->firstAtomInCell); 
}


PAIR_ITERATOR pfs_newIter(int iAtom, PAIR_FINDER_SIMPLE* this)
{
   PAIR_ITERATOR iter;
   
   iter.ii = iAtom;
   iter.jj = -1; // assure ii != jj
   iter.advance = pfs_advanceIter;
   PAIR_FINDER_SIMPLE_STATE* state = ddcMalloc(sizeof(PAIR_FINDER_SIMPLE_STATE));
   state->iCell = 0;
   state->nCells = pfs_mkCellList(iAtom, state->cellList, this);
   state->grid = this;
   iter.pairFinderState = state;
   iter.advance(&iter);
   return iter;
}



/** Three responsibilities:
 *  2.  Advance regular iterator
 *  3.  Clean up when we become an end iterator
 *      - iter->state needs to be freed
 *      - set iter->jj = iter->ii
 *
 */
void pfs_advanceIter(PAIR_ITERATOR* iter)
{
   PAIR_FINDER_SIMPLE_STATE* state = iter->pairFinderState;

   for(;state->iCell<state->nCells; ++state->iCell)
   {
      do
      {
	 if (iter->jj < 0)
	    iter->jj = state->grid->firstAtomInCell[state->cellList[state->iCell]];
	 else
	    iter->jj = state->grid->nextAtomInCell[iter->jj];
	 
	 if (iter->ii == iter->jj)
	    continue;
	 if (iter->jj >= 0)
	 {
	    double x = state->grid->r[iter->ii].x - state->grid->r[iter->jj].x;
	    double y = state->grid->r[iter->ii].y - state->grid->r[iter->jj].y;
	    double z = state->grid->r[iter->ii].z - state->grid->r[iter->jj].z;
	    double r2 = x*x+y*y+z*z;
	    if (r2 > state->grid->rCut2)
	    {
	       state->grid->nearestImage(&x, &y, &z);
	       r2 = x*x + y*y + z*z;
	    }
	    iter->rx = x;
	    iter->ry = y;
	    iter->rz = z;
	    iter->r2 = r2;
	    if (iter->r2 <= state->grid->rCut2)
	       return;
	 }
      }
      while (iter->jj >= 0);
   }

   // We have reached the end of the possible pairs.  Free resources to
   // avoid memory leaks and turn this iterator into an end iterator.
   iter->jj = iter->ii;
   ddcFree(state);
}

/**
 *  Make a list of all of the cells in the grid that particle iAtom
 *  could possibly interact with.  Since the grid spacing was chosen to
 *  be larger than rCut in most cases there will be 27 such cells.
 *  However, if the cutoff exceeds 1/3 of the box size along an axis
 *  there will be only 2 grid cells in that direction.  In that case the
 *  periodic reduction would cause us to add the same cell to the list
 *  twice.  Hence when adding a cell to the list we always scan the list
 *  to make sure that it isn't already there.
 */
int pfs_mkCellList(unsigned iAtom, int* cellList, PAIR_FINDER_SIMPLE* this)
{
   int ix = (this->r[iAtom].x - this->rMin.x) * this->drInv.x;
   int iy = (this->r[iAtom].y - this->rMin.y) * this->drInv.y;
   int iz = (this->r[iAtom].z - this->rMin.z) * this->drInv.z;
   int nx = this->nCell.x;
   int ny = this->nCell.y;
   int nz = this->nCell.z;

   while (ix<0) ix += nx;  while (ix>=nx) ix -= nx;
   while (iy<0) iy += ny;  while (iy>=ny) iy -= ny;
   while (iz<0) iz += nz;  while (iz>=nz) iz -= nz;
   int index = ix + nx*(iy + ny*iz);
   assert(index < nx*ny*nz);
   cellList[0] = index;

   int count = 1;
   for (int ii=-1; ii<2; ++ii)
      for (int jj=-1; jj<2; ++jj)
	 for (int kk=-1; kk<2; ++kk)
	 {
	    if (ii==0 && jj==0 &&kk==0) continue;

	    int x = ix + ii;
	    int y = iy + jj;
	    int z = iz + kk;
	    while (x<0) x += nx;  while (x>=nx) x -= nx;
	    while (y<0) y += ny;  while (y>=ny) y -= ny;
	    while (z<0) z += nz;  while (z>=nz) z -= nz;
	    int index = x + nx*(y + ny*z);
	    assert(index < nx*ny*nz);

	    cellList[count++] = index;

	    for (int tt=0; tt<count-1; ++tt)
	       if (index == cellList[tt])
	       {
		  --count;
		  break;
	       }
	 }
   return count;
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
