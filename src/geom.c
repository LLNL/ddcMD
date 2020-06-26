#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <limits.h>
#include "geom.h"
#include "box.h"
#include "error.h"
#include "ddcMalloc.h"
#include "heap.h"
#include "preduce.h"
#include "ptiming.h"

// To Do:
//
// The first step of the gridding algorithm is determining the spatial
// extent of the grid.  We can always do this by iterating the
// particles, doing a preduce to some center coordinate to make the
// particle coordinates compact, and computing min/max.  We can even
// store the compact coordinates so that we don't have to repeat the
// preduce when we later compute which grid cell the particle is in.
// However, this is a great deal of work that may be completely uneeded.
// We probably already have a good estimate for the extend of the grid
// because we know the radius of the domain that we are in.  Also, the
// coordinates may already be compact.  If both of these cases are true
// then we only need to stream once through the cooridinates to grid
// them.  We need no preduce and no storage for the particle
// coordinates.  This class should be written so that it can take
// advantage of those cases when possible.
//
// We could probably get substantial efficiency gains from detecting
// situations where the size of the grid has not changed and avoiding
// recomputing the offset arrays.
//
// In ddcFindPairs the PARTICLESET is build from ddc->particles.  It
// appears (from the fact that the preduce has been eliminated from the
// distance calculation at about line 149) that this is already a
// compact set.  Ensure that this is the case and then explore the
// resulting optimization posibilities.
//
// It would probably be faster to do the boxsep calculation and
// comparison in normalized coordinates.  However, that would require
// another set of Preduce functions that work on normalized coordinates.
// I don't feel like the hassle right now.  In fact, the whole finding
// of lists of neighbor boxes is a giant optimization target.  There are
// lots of possible improvements both in run time and storage space
// depending on the assumptions that are made.
//
//

static  enum GEOM_METHOD geomExpertSystem(GEOM* geom);
static void geomMethodDefault(GEOM* geom, PARTICLESET* ps);
static void geomMethodFast(GEOM* geom, PARTICLESET* ps);
static void geomMethodP7N1Ortho(GEOM* geom, PARTICLESET* ps);
static void geomNNBoxBox(GEOM *geom);
static void GeomDenseBox(GEOM*geom, PARTICLESET*ps);
static void geomBuildOffsetsFast(GEOM* geom);
void geomMethodP7N1Ortho(GEOM* geom, PARTICLESET* ps);
static THREE_VECTOR computeBoxSpan(
	THREE_MATRIX h, THREE_VECTOR* n0, THREE_VECTOR* n1, THREE_VECTOR* n2);

static const double _coarsenFactor = 1.2599; // (2.0)^(1/3)

/**
 *  Computes the shortest distance between b1 and b2.  This function is
 *  supposed to work for an arbitrary h matrix.
 *
 *  Returns distance squared.
 */
static double boxsep(const GEOMBOX *b1, const GEOMBOX *b2, const GEOM* geom)
{
	THREE_VECTOR rij;
	rij.x = b1->r.x - b2->r.x;
	rij.y = b1->r.y - b2->r.y;
	rij.z = b1->r.z - b2->r.z;
	backInBox_fast(&rij.x, &rij.y, &rij.z);

	THREE_VECTOR spanij;
	spanij.x = DOT(rij, geom->n0);
	spanij.y = DOT(rij, geom->n1);
	spanij.z = DOT(rij, geom->n2);

	rij.x = MAX(fabs(spanij.x) - geom->cellSpan.x, 0.0);
	rij.y = MAX(fabs(spanij.y) - geom->cellSpan.y, 0.0);
	rij.z = MAX(fabs(spanij.z) - geom->cellSpan.z, 0.0);
	double r2 = VSQ(rij);
	return r2; 
}


/** To do its work, GEOM needs a cutoff distance for pairs, the box
 * matrix, and its inverse.  We store pointers to the matrices so that
 * persistent GEOMs will automatically know about changes to the
 * h-matrix. */
GEOM* geom_init(double rCut, BOX_STRUCT *box)
{
	GEOM* this = ddcMalloc(sizeof(GEOM));
	this->method = GEOM_NONE;
   this->compBox = box; 
	this->h = &box->h0;
	this->hinv = &box->hinv;
   this->rCorner = box->reducedcorner; 
	this->boundary_condition = &box->pbc;
	// d, min, max, cellSpan, n0, n1, n2 uninitialized.
	this->minBoxSide = 0;
	this->rcut = rCut;
	this->small = 0;
	this->nx = 0;
	this->ny = 0;
	this->nz = 0;
	this->nbox = 0;
	this->box_nlist_update = 1;

	this->boxCapacity = 0;
	this->box = NULL;
	this->pinfoCapacity = 0;
	this->pinfo = NULL;
	this->r = NULL;
	this->boxOffsetCapacity = 0;
	this->boxOffset = NULL;
	return this;
}

void geom_free(GEOM* this)
{
	ddcFree(this->box);
	ddcFree(this->pinfo);
	ddcFree(this->boxOffset);
	ddcFree(this);
}

void geomMethodSet(GEOM* this, enum GEOM_METHOD method)
{
	this->method = method;
}



/** Do not use static variables or ExpandBuffers in this code.  We need
 * to be able to have separate, independent GEOM instances (no static
 * variables) and we need to be able to create a geom and throw it way
 * (no ExpandBuffers). */
static void geomNNBoxBox(GEOM *geom)
{
	const int maxBoxes = 27;

	int i,j;
//	static THREE_VECTOR **Rb, *R; 
	double r2; 

//	Rb = (THREE_VECTOR **) ExpandBuffers((void *)Rb, sizeof(THREE_VECTOR *), geom->nbox*27, 64, LOCATION("geomNNBoxBox"),"rb"); 
//	box_get(NULL,NEAREST_IMAGES,&R);
	
/*	if (!geom->box_nlist_update)  return;  */
	if (geom->nbox*maxBoxes > geom->boxOffsetCapacity)
	{
		geom->boxOffsetCapacity = 1.05*(geom->nbox)*maxBoxes;
		geom->boxOffset = (int*)
			ddcRealloc(geom->boxOffset, geom->boxOffsetCapacity * sizeof(int));
	}
	double rcut=geom->rcut;
	
//	for (i = 0; i < geom->nbox; i++) geom->box[i].R = Rb+maxBoxes*i;
	for (i = 0; i < geom->nbox; i++)
		geom->box[i].nlist = geom->boxOffset+maxBoxes*i;
	for (i = 0; i < geom->nbox; i++)
	{
		if (geom->box[i].first > -1) 
		{
//			k = (0+1)+ (0+1)*3 + (0+1)*9;
//			geom->box[i].R[geom->box[i].nn] = R+k;  
			geom->box[i].nlist[geom->box[i].nn++] = 0;

			for (j = 0; j < i; j++)
			{
				if (geom->box[j].first > -1) 
				{
					r2 = boxsep(geom->box+i,geom->box+j, geom);
					if (r2 < rcut*rcut) 
					{
//						preduceGetS(&S);  
//						ix = lrint (S.x); iy = lrint (S.y); iz = lrint (S.z);
//						k = (ix+1)+ (iy+1)*3 + (iz+1)*9;
						//geom->box[i].R[geom->box[i].nn] = R+k;   
//						k = (-ix+1)+ (-iy+1)*3 + (-iz+1)*9;
						//geom->box[j].R[geom->box[j].nn] = R+k;  
//						if (!geom->box_nlist_update)
						//{
							//int di,dj,ni,nj; 
							//ni = geom->box[i].nn; 
							//nj = geom->box[j].nn; 
							//di = geom->box[i].nlist[ni]; 
							//dj = geom->box[j].nlist[nj]; 
							//if (getRank(0) == 0 && ( di != (j-i) || dj != (i-j) )) printf("%d %d %d %d %d %d\n",ni,nj,di,dj,j-i,i-j); 
						//}
						geom->box[i].nlist[geom->box[i].nn++] = j -i;
						geom->box[j].nlist[geom->box[j].nn++] = i-j;
					}
				}
			}
		}
	}
	for (int ii=0; ii<geom->nbox; ++ii)
   {
		assert(geom->box[ii].nn <= maxBoxes);
    }
}

/** This routines handles building the offset arrays that give the
 *  neighbors of each box for one very specific, but very common case.
 *  This routine can be used when the following is true:
 *  - Each GEOM_BOX is at least rCut in size in its smallest dimension
 *  - We have a "dense" grid representation
 *  - The outermost "shell" of GEOM_BOXes contain only remote particles
 *    and therefore never need to access offsets.
 *
 *  To ensure that the latter condition is satisfied the offset pointers
 *  for all of the outer shell GEOM_BOXes will be set to NULL so that
 *  any attempt to use them will meet with and easily diagnosed failure.
 *
 *  In this case we know that we need to consider only the 26 nearest
 *  neighbor boxes.  We also know that we won't have to deal with any
 *  boundary issues.  Hence every GEOM_BOX that we will actually use has
 *  *exactly* the same offsets.  So, we can construct that array once
 *  and point all the GEOM_BOXes at that array.  Presto, instant time
 *  and memory savings.
 */
void geomBuildOffsetsFast(GEOM* geom)
{
	unsigned nx = geom->nx;
	unsigned ny = geom->ny;
	
	if (geom->boxOffsetCapacity < 27)
	{
		geom->boxOffsetCapacity = 27;
		geom->boxOffset = (int*) ddcMalloc(27*sizeof(int));
	}
	
	{
		unsigned cnt = 0;
		for (int ix=-1; ix<=1; ++ix)
			for (int iy=-1; iy<=1; ++iy)
				for (int iz=-1; iz<=1; ++iz)
				{
					geom->boxOffset[cnt++] = ix + nx*(iy + ny*iz);
				}
	}
	
	
	for (int ii=0; ii<geom->nbox; ++ii)
	{
		geom->box[ii].nn = 0;
		geom->box[ii].nlist = NULL;
	}

	for (int ix=1; ix<geom->nx-1; ++ix)
		for (int iy=1; iy<geom->ny-1; ++iy)
			for (int iz=1; iz<geom->nz-1; ++iz)
			{
				unsigned index = ix + nx*(iy + ny*iz);
				geom->box[index].nn = 27;
				geom->box[index].nlist = geom->boxOffset;
			}
}

void    imageParticle(GEOM *geom)
{
   int nx = geom->nx; 
   int ny = geom->ny; 
   int nz = geom->nz; 
   printf("%d %d %d\n",nx,ny,nz); 
   for (int i=0;i<nx*ny*nz;i++)
   {
      GEOMBOX *box = geom->box+i; 
      int type=box->type ;
      if (type != 13)
      {
      int ix = -(type%3-1); 
      int iy = -((type/3)%3-1); 
      int iz = -((type/9)%3-1); 
      int offset  = ix*(nx-2) + iy*nx*(ny-2) + iz*nx*ny*(nz-2); 
      GEOMBOX *imageBox = box+offset; 
      THREE_VECTOR *r = geom->r; 
      int j = imageBox->first;
      while (j != -1)
      {
         printf("%d %d %d %d %d %d %f %f %f\n", box->ix,box->iy,box->iz,imageBox->ix,imageBox->iy,imageBox->iz,r[j].x ,r[j].y, r[j].z); 
         j = geom->pinfo[j].next;
      }
      }
   }
   exit(0); 
}


/**
 *  First we preduce the particles with respect the "center" coordinate
 *  provided by the particleset to get a compact collection of
 *  coordindates.  In the same loop we also convert the coordinates to
 *  normalized box coordinates.
 *
 *  Second, we iterate the compacted particle positions to determine the
 *  grid parameters.  We find the min and max coord in each direction as
 *  well as the spacing.  Note that the grid spacing, d, is in normalized
 *  box units.
 */
void GeomBox(GEOM*geom, PARTICLESET*ps)
{
   double rcut = MAX(geom->rcut, geom->minBoxSide);
   assert(rcut > 0.0);
   assert(geom->h != NULL);
   assert(geom->hinv != NULL);

   int nParticles = *(ps->number_particles);
   unsigned rBlk;
   geom->r = (THREE_VECTOR*) heapGet(&rBlk);
   heapEndBlock(rBlk, nParticles * sizeof(THREE_VECTOR));

   if (nParticles > geom->pinfoCapacity)
   {
      geom->pinfoCapacity = 1.05 * nParticles + 512;
      geom->pinfo = (GEOM_PINFO*)
         ddcRealloc(geom->pinfo, sizeof(GEOM_PINFO)*geom->pinfoCapacity);
   }

   double* xptr = ps->xptr;
   double* yptr = ps->yptr;
   double* zptr = ps->zptr;
   unsigned stride = ps->stridedouble;

   for (int ii=0; ii<nParticles; ++ii)
   {
      double x = *(double*)(((char*)xptr) + ii*stride)-ps->center->x;
      double y = *(double*)(((char*)yptr) + ii*stride)-ps->center->y;
      double z = *(double*)(((char*)zptr) + ii*stride)-ps->center->z;
      backInBox_fast(&x, &y, &z);  
      VSET(geom->r[ii], x, y, z);
      geom->r[ii] = matrix_vector(*geom->hinv, geom->r[ii]);
   }

   THREE_VECTOR min = geom->r[0];
   THREE_VECTOR max = geom->r[0];
   for (int ii=1; ii<nParticles; ++ii)
   {
      THREE_VECTOR r = geom->r[ii];
      if (min.x > r.x) min.x = r.x;
      if (min.y > r.y) min.y = r.y;
      if (min.z > r.z) min.z = r.z;
      if (max.x < r.x) max.x = r.x;
      if (max.y < r.y) max.y = r.y;
      if (max.z < r.z) max.z = r.z;
   }
   geom->max = max; 
   geom->min = min; 

   geom->method = geomExpertSystem(geom);
   switch (geom->method)
   {
      case GEOM_DEFAULT:
         profile(BOX_DEFAULT, START);
         geomMethodDefault(geom, ps);
         profile(BOX_DEFAULT, END);
         break;
      case GEOM_FAST:
         profile(BOX_FAST, START);
         geomMethodFast(geom, ps);
         profile(BOX_FAST, END);
         break;
      case GEOM_P7N1ORTHO:
         geomMethodP7N1Ortho(geom, ps);
         geomBuildOffsetsFast(geom); 
         imageParticle(geom);
         break;
      default:
         assert(0 == 1);
   }
   heapFree(rBlk);
   geom->r = NULL;
}

/** Coordinate assigned to each cell is a real coordinate */
static void GeomDenseBox(GEOM*geom, PARTICLESET*ps)
{
   unsigned number_particles = *(ps->number_particles);
   unsigned number_local = *(ps->number_local);
   THREE_VECTOR min = geom->min; 
   int nx=geom->nx;
   int ny=geom->ny;
   int nz=geom->nz;
   THREE_VECTOR d=geom->d; 
   int nbox = nx*ny*nz;
   geom->nbox = nbox;
   if (nbox+1 > geom->boxCapacity)
   {
      geom->boxCapacity = 1.05*nbox+64;
      geom->box = ddcRealloc(geom->box, geom->boxCapacity * sizeof(GEOMBOX));
   }
   THREE_MATRIX h = *geom->h;
   for (int iz=0; iz<nz; iz++)
      for (int iy=0; iy<ny; iy++)
         for (int ix=0; ix<nx; ix++)
         {
            THREE_VECTOR r;
            r.x = ix*d.x + min.x  ;
            r.y = iy*d.y + min.y  ;
            r.z = iz*d.z + min.z  ;
            r = matrix_vector(h, r);
            int ll = ix + nx*(iy+ny*iz);
            int kk = 13; 
            if (ix == 0) kk-=1; if (ix == nx-1) kk+=1; 
            if (iy == 0) kk-=3; if (iy == ny-1) kk+=3; 
            if (iz == 0) kk-=9; if (iz == nz-1) kk+=9; 
            geom->box[ll].ix = ix; 
            geom->box[ll].iy = iy; 
            geom->box[ll].iz = iz; 
            geom->box[ll].key = ll ;
            geom->box[ll].index = ll ;
            geom->box[ll].type = kk ;
            geom->box[ll].r = r;
            geom->box[ll].first = -1;
            geom->box[ll].firstlocal = -1;
         }

   for (unsigned i = 0; i < number_particles; i++)
   {
      THREE_VECTOR r = geom->r[i]; 
      int ix = (int)((r.x - min.x)/d.x);
      int iy = (int)((r.y - min.y)/d.y);
      int iz = (int)((r.z - min.z)/d.z);
      ix = MAX(MIN(ix, nx - 1), 0);
      iy = MAX(MIN(iy, ny - 1), 0);
      iz = MAX(MIN(iz, nz - 1), 0);
      /* 		if ( (ix==0 || iy==0 || iz==0 || */
      /* 				ix==nx-1 || iy==ny-1 || iz ==nz-1) && */
      /* 			  i<number_local) */
      /* 			printf("Warning: local particle on grid edge\n"); */
      int bindex = ix + nx*(iy + ny*iz);

      geom->pinfo[i].box = geom->box + bindex;
      geom->pinfo[i].next = geom->box[bindex].first;
      geom->box[bindex].first = i;
      if (i < number_local)
         geom->box[bindex].firstlocal = i;
   }

   for (int i = 0; i < nbox; i++)
   {
      geom->box[i].nn = 0;
   }
}

PARTICLESET* ParticleSet(
      PARTICLESET* p, double *xptr, double *yptr, double *zptr, int stridedouble, gid_type *global_index,
      int stridelong64, int *type, int strideint, unsigned* number_particles, unsigned* number_local,
      const THREE_VECTOR* center)
{
   if (p == NULL) p = ddcMalloc(1*sizeof(PARTICLESET));
   p->xptr = xptr;
   p->yptr = yptr;
   p->zptr = zptr;
   p->global_index = global_index;
   p->type = type;
   p->number_particles = number_particles;
   p->number_local = number_local;
   p->stridedouble = stridedouble;
   p->strideint = strideint;
   p->stridelong64 = stridelong64;
   p->center = center;
   return p;
}


/** Computes the "span" of the box with box matrix h.  Here, we define
 *  the span as the distance between the planes on opposite sides of the
 *  box. span.x is the distance between planes connected by a0 and so
 *  on.*/
THREE_VECTOR computeBoxSpan(THREE_MATRIX h, THREE_VECTOR* n0, THREE_VECTOR* n1, THREE_VECTOR* n2)
{
   THREE_VECTOR a0 = {1, 0, 0};
   THREE_VECTOR a1 = {0, 1, 0};
   THREE_VECTOR a2 = {0, 0, 1};
   a0 = matrix_vector(h, a0);
   a1 = matrix_vector(h, a1);
   a2 = matrix_vector(h, a2);

   double a0_len = sqrt(dot1(a0, a0));
   double a1_len = sqrt(dot1(a1, a1));
   double a2_len = sqrt(dot1(a2, a2));

   *n0 = cross(&a1, &a2);
   *n1 = cross(&a2, &a0);
   *n2 = cross(&a0, &a1);
   VSCALE((*n0), 1.0/(a1_len * a2_len));
   VSCALE((*n1), 1.0/(a2_len * a0_len));
   VSCALE((*n2), 1.0/(a0_len * a1_len));

   THREE_VECTOR bodyDiag;
   bodyDiag.x = a0.x + a1.x + a2.x;
   bodyDiag.y = a0.y + a1.y + a2.y;
   bodyDiag.z = a0.z + a1.z + a2.z;

   THREE_VECTOR span;
   span.x = fabs(dot1(*n0, bodyDiag));
   span.y = fabs(dot1(*n1, bodyDiag));
   span.z = fabs(dot1(*n2, bodyDiag));
   return span;
}

enum GEOM_METHOD geomExpertSystem(GEOM* geom)
{


   THREE_VECTOR n0, n1, n2;
   THREE_VECTOR span = computeBoxSpan(*geom->h, &n0, &n1, &n2);

   THREE_VECTOR extent;
   extent.x = geom->max.x - geom->min.x;
   extent.y = geom->max.y - geom->min.y;
   extent.z = geom->max.z - geom->min.z;

   double rcut = MAX(geom->rcut, geom->minBoxSide);

   //BOX_STRUCT *box = geom->compBox; 
   //if((box->itype == ORTHORHOMBIC)  && (box->pbc == 7) && (getSize(0) == 1)) 
   //   return GEOM_P7N1ORTHO; 
   if ( (1-extent.x) * span.x  > rcut && (1-extent.y) * span.y  > rcut && (1-extent.z) * span.z  > rcut ) 
      return GEOM_FAST; 
   return GEOM_DEFAULT;
}



void geomMethodDefault(GEOM* geom, PARTICLESET* ps)
{
   assert(geom->method == GEOM_DEFAULT);

   double rcut = MAX(geom->rcut, geom->minBoxSide);
   THREE_VECTOR extent;
   extent.x = geom->max.x - geom->min.x;
   extent.y = geom->max.y - geom->min.y;
   extent.z = geom->max.z - geom->min.z;

   THREE_VECTOR boxSpan =
      computeBoxSpan(*geom->h, &geom->n0, &geom->n1, &geom->n2);
   // smallest spacing (in normalized coords) such that cell span > rcut.
   THREE_VECTOR d; 
   d.x = rcut/boxSpan.x;
   d.y = rcut/boxSpan.y;
   d.z = rcut/boxSpan.z;

   // Nx, ny, and nz are kept as doubles to avoid possible integer
   // overflow in the computation of nx*ny*nz.
   double nx = floor(extent.x/d.x); nx = MAX(1.0, nx);
   double ny = floor(extent.y/d.y); ny = MAX(1.0, ny);
   double nz = floor(extent.z/d.z); nz = MAX(1.0, nz);

   // Reduce grid density if it is too fine.  Ensure loop still exits if
   // number_particles == 0
   while (nx*ny*nz -1 > *(ps->number_particles))
   {
      nx = floor(nx / _coarsenFactor); nx = MAX(1, nx);
      ny = floor(ny / _coarsenFactor); ny = MAX(1, ny);
      nz = floor(nz / _coarsenFactor); nz = MAX(1, nz);
   }

   assert(nx < INT_MAX);
   assert(ny < INT_MAX);
   assert(nz < INT_MAX);

   // cell size in norm coords
   d.x = extent.x/nx; 
   d.y = extent.y/ny;
   d.z = extent.z/nz;

   // cell span in real length units.
   geom->cellSpan.x = extent.x*boxSpan.x/nx;
   geom->cellSpan.y = extent.y*boxSpan.y/ny;
   geom->cellSpan.z = extent.z*boxSpan.z/nz;


   // This block of code appears to be part of an incomplete attempt to
   // avoid regenerating the nlist when we know that the answer will be
   // the same lists we already have
   /*
      if (nx == geom->nx+1 &&  nx > 5 )nx --; 
      if (ny == geom->ny+1 &&  ny > 5 )ny --; 
      if (nz == geom->nz+1 &&  nz > 5) nz --; 
      geom->box_nlist_update=0; 
      if ( nx != geom->nx || ny != geom->ny || nz != geom->nz ) 
      {
      geom->box_nlist_update=1; 
      }
      */
   // Another incomplete attempt to improve the capabilities of GEOM
   /*	
      geom->small = 0;
      if (nx < 3 || ny < 3 || nz < 3) geom->small = 1;
      */

   geom->nx = nx ;
   geom->ny = ny ;
   geom->nz = nz ;
   geom->d = d;
   GeomDenseBox(geom, ps);
   /* 	double rcut = geom->rcut; */
   /* 	double minspan = box_get_minspan(NULL);  */
   /* 	double boxspan =  nx*d.x;  */
   /* 	boxspan = MAX(boxspan,ny*d.y); */
   /* 	boxspan = MAX(boxspan,nz*d.z); */
   /* 	minspan -= boxspan;  */
   /*	if (minspan < rcut) geomNNRegPeriodic(geom);  */
   /* 	if (minspan < rcut || 1==1) geomNNBoxBox(geom);  */
   //	else geomNNRegFree(geom);   
   geomNNBoxBox(geom);
}

void geomMethodFast(GEOM* geom, PARTICLESET* ps)
{
   assert(geom->method == GEOM_FAST);

   double rcut = MAX(geom->rcut, geom->minBoxSide);

   THREE_VECTOR boxSpan =
      computeBoxSpan(*geom->h, &geom->n0, &geom->n1, &geom->n2);
   // smallest spacing (in normalized coords) such that cell span > rcut.
   geom->d.x = rcut/boxSpan.x;
   geom->d.y = rcut/boxSpan.y;
   geom->d.z = rcut/boxSpan.z;

   THREE_VECTOR extent;
   extent.x = geom->max.x - geom->min.x;
   extent.y = geom->max.y - geom->min.y;
   extent.z = geom->max.z - geom->min.z;

   // In corner cases such as no particles, or particles that all lie in
   // an axis normal plane an extent can be zero.  This will lead to
   // nastiness since we always want at least three boxes in each
   // direction.  So, ensure that the extent is at least the computed
   // spacing.  Then extent/geom->d >= 1.0.
   extent.x = MAX(extent.x, geom->d.x);
   extent.y = MAX(extent.y, geom->d.y);
   extent.z = MAX(extent.z, geom->d.z);

   double nx = ceil(extent.x/geom->d.x)+2;
   double ny = ceil(extent.y/geom->d.y)+2;
   double nz = ceil(extent.z/geom->d.z)+2;

   // cell span in real length units.
   geom->cellSpan.x = rcut;
   geom->cellSpan.y = rcut;
   geom->cellSpan.z = rcut;

   // Reduce grid density if it is too fine.  Ensure loop still exits if
   // number_particles == 0.  Nx, ny, and nz are doubles to prevent
   // potential integer overflow of product.
   while (nx*ny*nz - 27 > *(ps->number_particles))
   {
      geom->d.x *= _coarsenFactor; geom->cellSpan.x *= _coarsenFactor;
      geom->d.y *= _coarsenFactor; geom->cellSpan.y *= _coarsenFactor;
      geom->d.z *= _coarsenFactor; geom->cellSpan.z *= _coarsenFactor;

      nx = ceil(extent.x/geom->d.x)+2;
      ny = ceil(extent.y/geom->d.y)+2;
      nz = ceil(extent.z/geom->d.z)+2;
   }

   assert(nx > 2 && nx < INT_MAX);
   assert(ny > 2 && ny < INT_MAX);
   assert(nz > 2 && nz < INT_MAX);

   geom->min.x = (geom->min.x+geom->max.x - nx*geom->d.x)/2.0;
   geom->min.y = (geom->min.y+geom->max.y - ny*geom->d.y)/2.0;
   geom->min.z = (geom->min.z+geom->max.z - nz*geom->d.z)/2.0;

   geom->max.x = geom->min.x + nx*geom->d.x;
   geom->max.y = geom->min.y + ny*geom->d.y;
   geom->max.z = geom->min.z + nz*geom->d.z;

   geom->nx = nx ;
   geom->ny = ny ;
   geom->nz = nz ;
   GeomDenseBox(geom, ps);
   geomBuildOffsetsFast(geom);
}
void geomMethodP7N1Ortho(GEOM* geom, PARTICLESET* ps)
{
   double rcut = MAX(geom->rcut, geom->minBoxSide);
   THREE_VECTOR extent;
   extent.x = 1.0;
   extent.y = 1.0;
   extent.z = 1.0;

   //		computeBoxSpan(*geom->h, &geom->n0, &geom->n1, &geom->n2);
   // smallest spacing (in normalized coords) such that cell span > rcut.
   THREE_VECTOR d; 
   d.x = rcut/geom->h->xx;
   d.y = rcut/geom->h->yy;
   d.z = rcut/geom->h->zz;

   // Nx, ny, and nz are kept as doubles to avoid possible integer
   // overflow in the computation of nx*ny*nz.
   double nx = floor(extent.x/d.x); nx = MAX(1.0, nx);
   double ny = floor(extent.y/d.y); ny = MAX(1.0, ny);
   double nz = floor(extent.z/d.z); nz = MAX(1.0, nz);

   assert(nx < INT_MAX);
   assert(ny < INT_MAX);
   assert(nz < INT_MAX);

   // cell size in norm coords
   d.x = extent.x/nx; 
   d.y = extent.y/ny;
   d.z = extent.z/nz;

   extent.x += 2*d.x ; 
   extent.y += 2*d.y ; 
   extent.z += 2*d.z ; 
   nx+=2; 
   ny+=2; 
   nz+=2; 
   geom->min.x = geom->rCorner.x-d.x; 
   geom->min.y = geom->rCorner.y-d.y; 
   geom->min.z = geom->rCorner.y-d.z; 

   geom->max.x = geom->min.x + extent.x; 
   geom->max.y = geom->min.y + extent.y; 
   geom->max.z = geom->min.z + extent.z; 

   // cell span in real length units.
   geom->cellSpan.x = d.x*geom->h->xx;
   geom->cellSpan.y = d.y*geom->h->yy;
   geom->cellSpan.z = d.z*geom->h->zz;


   geom->nx = nx ;
   geom->ny = ny ;
   geom->nz = nz ;
   geom->d = d;
   GeomDenseBox(geom, ps);
   /* 	double rcut = geom->rcut; */
   /* 	double minspan = box_get_minspan(NULL);  */
   /* 	double boxspan =  nx*d.x;  */
   /* 	boxspan = MAX(boxspan,ny*d.y); */
   /* 	boxspan = MAX(boxspan,nz*d.z); */
   /* 	minspan -= boxspan;  */
   /*	if (minspan < rcut) geomNNRegPeriodic(geom);  */
   /* 	if (minspan < rcut || 1==1) geomNNBoxBox(geom);  */
   //	else geomNNRegFree(geom);   
   //	geomNNBoxBox(geom);
}



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//   
//  Everything below this point is dinosaur tracks.
//  All of this code was unused and unreachable when
//  I did the refactoring of the GEOM code to decouple
//  it from the NBR list in Sep 09.  It is preserved
//  here for possible future study and analysis since
//  some of the ideas are likely good.
//
//  -dfr Sep 09
//
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


#if 0

THREE_VECTOR nlattice[]={ { 1, 1, 1} ,{ 1,0, 1},{ 1,-1, 1},{0, 1, 1},{0,0, 1},{0,-1, 1},{-1, 1, 1},{-1,0, 1},{-1,-1, 1},{ 1, 1,0},{0, 1,0},{-1, 1,0},{ 1,0,0}}; 

typedef struct dn_str{ int x,y,z;double r;} DN; 
struct { int x,y,z;} dn[]={ { 1, 1, 1} ,{ 1,0, 1},{ 1,-1, 1},{0, 1, 1},{0,0, 1},{0,-1, 1},{-1, 1, 1},{-1,0, 1},{-1,-1, 1},{ 1, 1,0},{0, 1,0},{-1, 1,0},{ 1,0,0}}; 


void GeomSparseBox(GEOM*geom, PARTICLESET*ps, NPARTICLE*particles);
void geomNNSmallBox(GEOM *geom);

int boxsort(NPARTICLE *x, NPARTICLE *y)
{
   if (x->bindex < y->bindex) return -1; 
   if (x->bindex ==  y->bindex) return 0; 
   return 1;
}
int psort(NPARTICLE *x, NPARTICLE *y)
{
   if (x->pindex < y->pindex) return -1; 
   if (x->pindex ==  y->pindex) return 0; 
   return 1;
}
int  BoxKeyCompare(int *key, GEOMBOX *b)
{
   if (*key < b->key) return -1; 
   if (*key > b->key) return 1; 
   return 0; 
}

void GeomBoxNew(GEOM*geom, PARTICLESET*ps, NPARTICLE*particles)
{
   int i, j, l, ix, iy, iz, nx, ny, nz, nbox;
   int number_particles, number_local;
   int last,first; 
   THREE_VECTOR r, r0, min, max, d;
   double rcut;
   GEOMBOX *box;

   rcut = geom->rcut[0].value+(geom->deltaR);
   number_particles = *(ps->number_particles);
   number_local = *(ps->number_local);
   i = 0;
   r0 = particles[0].r0; 
   min = max = r0;
   for (i = 1; i < number_particles; i++)
   {
      r = particles[i].r0; 
      if (min.x > r.x) min.x = r.x;
      if (min.y > r.y) min.y = r.y;
      if (min.z > r.z) min.z = r.z;
      if (max.x < r.x) max.x = r.x;
      if (max.y < r.y) max.y = r.y;
      if (max.z < r.z) max.z = r.z;
   }
   d.x = max.x - min.x;
   d.y = max.y - min.y;
   d.z = max.z - min.z;
   if (d.x < rcut) d.x = 1.0001*rcut;
   if (d.y < rcut) d.y = 1.0001*rcut;
   if (d.z < rcut) d.z = 1.0001*rcut;
   geom->nx = nx = MAX((int)(d.x/(rcut)),1);
   geom->ny = ny = MAX((int)(d.y/(rcut)),1);
   geom->nz = nz = MAX((int)(d.z/(rcut)),1);
   geom->small = 0;
   if (nx < 3 || ny < 3 || nz < 3) geom->small = 1;
   d.x /= nx;
   d.y /= ny;
   d.z /= nz;
   geom->d=d; 
   geom->nbox = 0;
   l = 0;
   for (i = 0; i < number_particles; i++)
   {
      r = particles[i].r0; 
      ix = (int)((r.x - min.x)/d.x);
      iy = (int)((r.y - min.y)/d.y);
      iz = (int)((r.z - min.z)/d.z);
      ix = MAX(MIN(ix, nx - 1), 0);
      iy = MAX(MIN(iy, ny - 1), 0);
      iz = MAX(MIN(iz, nz - 1), 0);
      l = ix + nx*(iy + ny*iz);
      particles[i].pindex = i; 
      particles[i].bindex = l; 
   }
   qsort(particles,number_particles,sizeof(NPARTICLE),(int(*)(const void*,const void*))boxsort); 
   last =particles[0].bindex; 
   first =0 ;
   nbox = 0; 
   particles[0].bindex=nbox;
   for (i = 1; i < number_particles; i++)
   {
      if (particles[i].bindex != last )
      {
         geom->box = (GEOMBOX *) ExpandBuffers((void *)geom->box, sizeof(GEOMBOX), nbox+1, 64, LOCATION("geomBox"),"geom->box"); 
         ix = last % nx ; last /=  nx ; iy = last % ny; last /=  ny ; iz = last; 
         r.x = ix*d.x + min.x  ;
         r.y = iy*d.y + min.y  ;
         r.z = iz*d.z + min.z  ;
         geom->box[nbox].ix = ix; 
         geom->box[nbox].iy = iy; 
         geom->box[nbox].iz = iz; 
         geom->box[nbox].key = ix + nx*(iy+ny*iz);
         geom->box[nbox].index = nbox; 
         geom->box[nbox].r = r;
         last = particles[i].bindex; 
         first=i;
         nbox ++; 
      }
      particles[i].bindex=nbox;
   }
   geom->box = (GEOMBOX *) ExpandBuffers((void *)geom->box, sizeof(GEOMBOX), nbox+1, 64, LOCATION("geomBox"),"geom->box"); 
   ix = last % nx ; last /=  nx ; iy = last % ny; last /=  ny ; iz = last; 
   r.x = ix*d.x + min.x  ;
   r.y = iy*d.y + min.y  ;
   r.z = iz*d.z + min.z  ;
   geom->box[nbox].key = ix + nx*(iy+ny*iz);
   geom->box[nbox].index = nbox; 
   geom->box[nbox].r = r;
   nbox ++; 
   geom->nbox = nbox; 
   qsort(particles,number_particles,sizeof(NPARTICLE),(int(*)(const void*,const void*))psort); 
   for (i = 0; i < number_particles; i++)
   {
      particles[i].box = geom->box + particles[i].bindex ;
      particles[i].next =-1; 
   }
   for (i = 0; i < nbox; i++)
   {
      geom->box[i].nn = 0;
      geom->box[i].first = -1;
      geom->box[i].firstlocal = -1;
   }
   for (i = 0; i < number_particles; i++)
   {
      box = particles[i].box;
      particles[i].next = box->first;
      j = particles[i].pindex; 
      box->first = j;
      if (j < number_local) box->firstlocal = j;
   }
   geomNNBoxBox(geom);
}

/*
   void geomNNSmallBox(GEOM *geom)
   {
   int i,j,k,jmin;
   THREE_VECTOR d,R; 
   double rcut,r2; 
   unsigned int zstart[1024],iz;
   rcut=geom->rcut[0]; 
   d = geom->d; 
   for (i=0;i<geom->nz;i++) zstart[i] = -1; 
   for (i = 0; i < geom->nbox; i++) 
   {
   iz = geom->box[i].iz; 
   if (zstart[iz]==-1) zstart[iz] = i; 
   }
   for (i = 0; i < geom->nbox; i++)
   {
   if (geom->box[i].first > -1) 
   {
   geom->box[i].nlist[geom->box[i].nn++] = 0;

   iz = geom->box[i].iz; 
   jmin=0; 
   if (iz > 0) jmin=zstart[iz-1];
   if (jmin == -1) jmin = zstart[iz];
   jmin=0; 
   for (j = jmin; j < i; j++)
   {
   if (geom->box[j].first > -1) 
   {
   r2 = boxsep(geom->box+i,geom->box+j,&d);
   if (r2 < rcut*rcut) 
   {
   preduceGetR(&R);  
   geom->box[i].R[geom->box[i].nn] = R;  
   geom->box[i].nlist[geom->box[i].nn++] = j -i;
   R.x *= -1.0;R.y *= -1.0; R.z *= -1.0; 
   geom->box[j].R[geom->box[j].nn] = R; 
   geom->box[j].nlist[geom->box[j].nn++] = i-j;
   }
   }
   }
#if xxxxx
if (iz==(geom->nz-1) &&  iz > 0)
{
j=0; 
while (geom->box[j].iz ==0) 
{
if (geom->box[j].first > -1) 
{
r2 = boxsep(geom->box+i,geom->box+j,&d);
if (r2 < rcut*rcut) 
{
preduceGetS(&R);  
ix = lrint (R.x);
iy = lrint (R.y);
iz = lrint (R.z);
k = (ix+1)+ (iy+1)*3 + (iz+1)*9;
geom->box[i].R[geom->box[i].nn] = R;  
geom->box[i].nlist[geom->box[i].nn++] = j-i;
R.x *= -1.0;R.y *= -1.0; R.z *= -1.0; 
geom->box[j].R[geom->box[j].nn] = R; 
geom->box[j].nlist[geom->box[j].nn++] = i-j;
}
}
j++;
}
}
#endif
}
}
}
*/

void GeomSparseBox(GEOM*geom, PARTICLESET*ps, NPARTICLE*particles)
{
   int i, j, l, ix, iy, iz, nx, ny, nz, nbox;
   int number_particles, number_local;
   int last,first; 
   THREE_VECTOR r, min, max, d;
   double rcut;
   GEOMBOX *box;

   rcut = geom->rcut[0].value+(geom->deltaR);
   number_particles = *(ps->number_particles);
   number_local = *(ps->number_local);
   min = geom->min; 
   max = geom->max; 
   nx=geom->nx;
   ny=geom->ny;
   nz=geom->nz;
   d=geom->d; 
   geom->nbox = 0;
   l = 0;
   for (i = 0; i < number_particles; i++)
   {
      r = particles[i].r0; 
      ix = (int)((r.x - min.x)/d.x);
      iy = (int)((r.y - min.y)/d.y);
      iz = (int)((r.z - min.z)/d.z);
      ix = MAX(MIN(ix, nx - 1), 0);
      iy = MAX(MIN(iy, ny - 1), 0);
      iz = MAX(MIN(iz, nz - 1), 0);
      l = ix + nx*(iy + ny*iz);
      particles[i].pindex = i; 
      particles[i].bindex = l; 
   }
   qsort(particles,number_particles,sizeof(NPARTICLE),(int(*)(const void*,const void*))boxsort); 
   last =particles[0].bindex; 
   first =0 ;
   nbox = 0; 
   particles[0].bindex=nbox;
   for (i = 1; i < number_particles; i++)
   {
      if (particles[i].bindex != last )
      {
         geom->box = (GEOMBOX *) ExpandBuffers((void *)geom->box, sizeof(GEOMBOX), nbox+1, 64, LOCATION("geomBox"),"geom->box"); 
         ix = last % nx ; last /=  nx ; iy = last % ny; last /=  ny ; iz = last; 
         r.x = ix*d.x + min.x  ;
         r.y = iy*d.y + min.y  ;
         r.z = iz*d.z + min.z  ;
         geom->box[nbox].ix = ix; 
         geom->box[nbox].iy = iy; 
         geom->box[nbox].iz = iz; 
         geom->box[nbox].key = ix + nx*(iy+ny*iz);
         geom->box[nbox].index = nbox; 
         geom->box[nbox].r = r;
         last = particles[i].bindex; 
         first=i;
         nbox ++; 
      }
      particles[i].bindex=nbox;
   }
   geom->box = (GEOMBOX *) ExpandBuffers((void *)geom->box, sizeof(GEOMBOX), nbox+1, 64, LOCATION("geomBox"),"geom->box"); 
   ix = last % nx ; last /=  nx ; iy = last % ny; last /=  ny ; iz = last; 
   r.x = ix*d.x + min.x  ;
   r.y = iy*d.y + min.y  ;
   r.z = iz*d.z + min.z  ;
   geom->box[nbox].key = ix + nx*(iy+ny*iz);
   geom->box[nbox].index = nbox; 
   geom->box[nbox].r = r;
   nbox ++; 
   geom->nbox = nbox; 
   qsort(particles,number_particles,sizeof(NPARTICLE),(int(*)(const void*,const void*))psort); 
   for (i = 0; i < number_particles; i++)
   {
      particles[i].box = geom->box + particles[i].bindex ;
      particles[i].next =-1; 
   }
   for (i = 0; i < nbox; i++)
   {
      geom->box[i].nn = 0;
      geom->box[i].first = -1;
      geom->box[i].firstlocal = -1;
   }
   for (i = 0; i < number_particles; i++)
   {
      box = particles[i].box;
      particles[i].next = box->first;
      j = particles[i].pindex; 
      box->first = j;
      if (j < number_local) box->firstlocal = j;
   }
   geomNNBoxBox(geom);
}

/*
   void geomNNFreeBox(GEOM *geom)
   {
   unsigned int key,jkey; 
   int jx,jy,jz,ix,iy,iz; 
   int i,j,l,nx,ny,nz;
   GEOMBOX *jbox; 
   THREE_VECTOR d,R; 
   double rcut,r2; 
   nx = geom->nx;
   ny = geom->ny;
   nz = geom->nz;
   rcut=geom->rcut[0]; 
   d = geom->d; 
   for (i = 0; i < geom->nbox; i++)
   {
   if (geom->box[i].first > -1) 
   {
   key = geom->box[i].key; 
   ix = geom->box[i].ix;  
   iy = geom->box[i].iy;  
   iz = geom->box[i].iz;  
   geom->box[i].nlist[geom->box[i].nn++] = geom->box + i;
   for (l = 0; l < 13; l++)
   {
   jx = (ix + dn[l].x+nx)%nx ;
   jy = (iy + dn[l].y+ny)%ny;
   jz = (iz + dn[l].z+nz)%nz ;
   if (!(jx == -1 || jx ==nx || jy==-1 || jy ==ny || jz==-1 || jz == nz))
   {
   jkey = jx + nx*(jy+ny*jz); 
   jbox = bsearch(&jkey,geom->box,geom->nbox,sizeof(GEOMBOX),(int(*)(const void*,const void*))BoxKeyCompare); 
   if (jbox != NULL) 
   {
   j = jbox->index; 
   if (geom->box[j].first > -1) 
   {

   r2 = boxsep(geom->box+i,geom->box+j,&d);
   if (r2 < rcut*rcut) 
   {
   preduceGetR(&R);  
   geom->box[i].R[geom->box[i].nn] = R;  
   R.x *= -1.0;R.y *= -1.0; R.z *= -1.0; 
   geom->box[j].R[geom->box[j].nn] = R; 
   geom->box[i].nlist[geom->box[i].nn++] = geom->box + j;
   geom->box[j].nlist[geom->box[j].nn++] = geom->box + i;
   }
   }
   }
   }
   }
   }
   }
   }
   */

int sortdn(DN *a , DN *b)
{
   double eps = 1e-8;
   if ( a->r < b->r-eps ) return -1 ;
   if ( a->r > b->r+eps ) return 1 ;
   return 0; 
}

void geomNNRegFree(GEOM *geom)
{
#define L 1 
#define LL   (2*L+1)*(2*L+1)*(2L+1)
   int i,j,k,l,n,ix,iy,iz,jx,jy,jz,nx,ny,nz,type;
   double x,y,z,r,dx,dy,dz;
   DN dn[LL]; 
   static THREE_VECTOR *R[LL] , RR={0.0,0.0,0.0}; 
   static int nn[LL],list[LL][LL]; 
   nx = geom->nx; 
   ny = geom->ny; 
   nz = geom->nz; 
   dx = geom->d.x; 
   dy = geom->d.y; 
   dz = geom->d.z; 

   k=0; 
   for (jz=-L;jz<=L;jz++)
      for (jy=-L;jy<=L;jy++)
         for (jx=-L;jx<=L;jx++)
         {
            if ((jx*jx + jy*jy + jz*jz) > 0) 
            {
               x = MAX(fabs(jx*dx) - dx, 0.0);
               y = MAX(fabs(jy*dy) - dy, 0.0);
               z = MAX(fabs(jz*dz) - dz, 0.0);
               r = sqrt(x*x + y*y + z*z);
               dn[k].x = jx;
               dn[k].y = jy;
               dn[k].z = jz;
               dn[k].r = r;
               k++; 
            }

         }
   qsort(dn, k, sizeof(DN), (int(*)(const void*,const void*))sortdn);
   l=0; 
   for (iz=-L;iz<=L;iz++)
      for (iy=-L;iy<=L;iy++)
         for (ix=-L;ix<=L;ix++)
         {
            n=0; 
            list[l][n++]=0; 
            for (i=0;i<k;i++)
            {
               jx = dn[i].x;
               jy = dn[i].y;
               jz = dn[i].z;
               if (ix+jx > -2 && ix+jx < 2)
                  if (iy+jy > -2 && iy+jy < 2)
                     if (iz+jz > -2 && iz+jz < 2)
                     {
                        j =  jx + nx*(jy + ny*jz) ;
                        list[l][n++]=j; 
                     }
            }
            R[l]=&RR; 
            nn[l++]=n; 
         }
   /*
      for (iz=-L;iz<=L;iz++)
      for (iy=-L;iy<=L;iy++)
      for (ix=-L;ix<=L;ix++)
      {
      n=0; 
      list[l][n++]=0; 
      for (jz=-L;jz<=L;jz++)
      for (jy=-L;jy<=L;jy++)
      for (jx=-L;jx<=L;jx++)
      {
      if (ix+jx >= -L && ix+jx <= L)
      if (iy+jy >= -L && iy+jy <= L)
      if (iz+jz >= -L && iz+jz <= L)
      if ((jx*jx + jy*jy + jz*jz) > 0) 
      {
      j =  jx + nx*(jy + ny*jz) ;
      list[l][n++]=j; 
      }
      }
      for (i=0;i<13;i++)
      {


      jx = dn[i].x;
      jy = dn[i].y;
      jz = dn[i].z;
      if (ix+jx > -2 && ix+jx < 2)
      if (iy+jy > -2 && iy+jy < 2)
      if (iz+jz > -2 && iz+jz < 2)
      {
      j =  jx + nx*(jy + ny*jz) ;
      list[l][n++]=j; 
      }
      jx = -dn[i].x;
      jy = -dn[i].y;
      jz = -dn[i].z;
      if (ix+jx > -2 && ix+jx < 2)
      if (iy+jy > -2 && iy+jy < 2)
      if (iz+jz > -2 && iz+jz < 2)
      {
      j =  jx + nx*(jy + ny*jz) ;
      list[l][n++]=j; 
      }
      }
      R[l]=&RR; 
      nn[l++]=n; 
      }
      */
   for (i = 0; i < geom->nbox; i++)
   {
      type = geom->box[i].type;
      geom->box[i].R = R;  
      geom->box[i].nlist= list[type];
      geom->box[i].nn= nn[type];
   }
}
void geomNNRegPeriodic(GEOM *geom)
{
   int i,j,k,l,n,ix,iy,iz,jx,jy,jz,kx,ky,kz,nx,ny,nz,type;
   static THREE_VECTOR *R[27][27] ,*RR; 
   static int nn[27],list[27][27]; 
   nx = geom->nx; 
   ny = geom->ny; 
   nz = geom->nz; 
   l=0; 
   box_get(NULL,NEAREST_IMAGES,&RR);
   for (iz=-1;iz<=1;iz++)
      for (iy=-1;iy<=1;iy++)
         for (ix=-1;ix<=1;ix++)
         {
            n=0; 
            k=13; 
            R[l][n]= RR+k;
            list[l][n++]=0; 
            for (i=0;i<13;i++)
            {
               jx = dn[i].x;
               jy = dn[i].y;
               jz = dn[i].z;
               kx=ky=kz=0; 
               if (ix+jx == -2 ) kx=-(jx = 1); 
               if (ix+jx ==  2 ) kx=-(jx =-1) ; 
               if (iy+jy == -2 ) ky=-(jy = 1) ; 
               if (iy+jy ==  2 ) ky=-(jy =-1) ; 
               if (iz+jz == -2 ) kz=-(jz = 1) ; 
               if (iz+jz ==  2 ) kz=-(jz =-1) ; 
               j =  jx + nx*(jy + ny*jz) ;
               k = kx + 3*ky + 9*kz; 
               R[l][n]= RR+k;
               list[l][n++]=j; 
               jx = -dn[i].x;
               jy = -dn[i].y;
               jz = -dn[i].z;
               kx=ky=kz=0; 
               if (ix+jx == -2 ) kx=-(jx = 1); 
               if (ix+jx ==  2 ) kx=-(jx =-1) ; 
               if (iy+jy == -2 ) ky=-(jy = 1) ; 
               if (iy+jy ==  2 ) ky=-(jy =-1) ; 
               if (iz+jz == -2 ) kz=-(jz = 1) ; 
               if (iz+jz ==  2 ) kz=-(jz =-1) ; 
               j =  jx + nx*(jy + ny*jz) ;
               k = kx + 3*ky + 9 *kz; 
               R[l][n]= RR+k;
               list[l][n++]=j; 
            }
            nn[l++]=n; 
         }
   for (i = 0; i < geom->nbox; i++)
   {
      type = geom->box[i].type;
      geom->box[i].R = R[type] ;  
      geom->box[i].nlist= list[type];
      geom->box[i].nn= nn[type];
   }
}

#endif

/* Local Variables: */
/* tab-width: 3 */
/* End: */
