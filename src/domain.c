#include "domain.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>

#include <unistd.h>
#include <limits.h>

#include "object.h"
#include "ddc.h"
#include "box.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include "utilities.h"
#include "preduce.h"
#include "mpiTypes.h"

#include "state.h"
#include "system.h"

#include "rectimer.h"
#include "functimer.h"


static void cubic(DDC* ddc);
static void bcc  (DDC* ddc);
static void fcc  (DDC* ddc);
static void hexagonal(DDC*ddc);
static void nonlattice(DDC* ddc);

static void recursive_bisection_domset(DDC *ddc);

static void generateLattice(DDC* ddc, int nBasis, THREE_VECTOR* basis);
static void mkBinWidth(int nr, double* dr, double* binWidth, int periodic);

void domainset_init(DOMAINSET* dset)
{
   int ndomains;
   MPI_Comm_size( COMM_LOCAL, &ndomains );
   dset->size     = ndomains;
   dset->domains  = (DOMAINX *) ddcMalloc(dset->size*sizeof(DOMAINX));

   int id;
   MPI_Comm_rank(COMM_LOCAL, &id);
   dset->local_id=id;
   dset->local_domain  =  dset->domains+id;
}

//  From my current analysis of this function, if we don't care about
//  setting the domain radii and bounding boxes (and I don't think we
//  do) then all of the nonlattice/qhull code *except* the call to
//  computeLocalVoronoiInfo can go away.  And that call should probably
//  go elsewehere.  In particular, we can get rid of the call to
//  computeLocalVoronoiRadius

void domain_set(DDC*ddc)
{
   profile(DOMAIN_SET, START);
   assert( ddc->domains.domains!=NULL );
   assert( ddc->domains.size==ddc->size );

   switch (ddc->lattice)
   {
      case CUBIC:
         cubic(ddc);
         break;
      case BCC:
         bcc(ddc);
         break;
      case FCC:
         fcc(ddc);
         break;
      case HEXAGONAL:
         hexagonal(ddc);
         break;
      case NONLATTICE:
         nonlattice(ddc);
         break;
      case RECURSIVE_BISECTION:
         /*printf("Choosing recursive bisection branch in domain_set() switch.\n");*/
         recursive_bisection_domset(ddc);
         break;
      case NONE:
         if (getRank(0) ==0)
            printf("Number of task must be in the set {lx*lx*lz,2*lx*ly*lz^3,4*lx*ly*lz lx=%d ly=%d lz=%d\n",
                  ddc->lx,ddc->ly,ddc->lz);
         abortAll(1);
         break;
      default:
         assert(1==0);
         break;
   }
   profile(DOMAIN_SET, END);
}

int domain_possibleRemote(const DOMAINX d, THREE_VECTOR v, double rcut)
{
   v.x -= d.center.x;
   v.y -= d.center.y;
   v.z -= d.center.z;
   nearestImage(&v.x, &v.y, &v.z);
   switch(d.shape)
   {
      double r2,R;
      case SPHERE:

      r2 = v.x*v.x + v.y*v.y + v.z*v.z;
      R = d.radius + rcut;
      if (r2 < (R*R)) return 1;
      break;
      case BRICK:
      if ( (fabs(v.x) <= d.extent.x + rcut) &&
            (fabs(v.y) <= d.extent.y + rcut) &&
            (fabs(v.z) <= d.extent.z + rcut) )
         return 1;
      break;
   }
   return 0;
}
/**
 *  Populates u with two vectors.
 *  u[0] is a unit vector in the direction of c1 to c2.
 *  u[1] is the midpoint between c1 and c2.
 *  Returns the distance from c1 to c2
 */
double domainset_separatingPlane(DOMAINX* local_domain, DOMAINX* remote_domain, THREE_VECTOR *u)
{
   const THREE_VECTOR *c1 = &(local_domain->center);
   const THREE_VECTOR* c2 = &(remote_domain->center);
   double ui,dist; 
   u[0].x = c2->x - c1->x; 
   u[0].y = c2->y - c1->y; 
   u[0].z = c2->z - c1->z; 
   backInBox_fast(&u[0].x, &u[0].y, &u[0].z);
   u[1].x = 0.5*u[0].x +c1->x; 
   u[1].y = 0.5*u[0].y +c1->y; 
   u[1].z = 0.5*u[0].z +c1->z; 
   dist = sqrt(VSQ(u[0])); 
   ui = 1.0/dist; 
   u[0].x *= ui; 
   u[0].y *= ui; 
   u[0].z *= ui; 
   return  dist; 
}
int domain_possibleRemote_by_plane(DDC*ddc, THREE_VECTOR *u, THREE_VECTOR v, double *dist)
{
   double d;
   v.x -= u[1].x; 
   v.y -= u[1].y; 
   v.z -= u[1].z; 
   backInBox_fast(&v.x, &v.y, &v.z);
   *dist = d =-(u[0].x*v.x + u[0].y*v.y + u[0].z*v.z);
   if (d < ddc->rcut) return 1; 
   return 0;
}


// determine to which domain a particle at position v belong to
// using the shortest distance to domain centers
int domainset_particle(DOMAINSET* this, THREE_VECTOR v, int n, int *list)
{
   backInBox_fast(&v.x, &v.y, &v.z);

   int jmin=-1;
   double r2min = 1e307;

   for (int j = 0; j < n; j++)
   {
      THREE_VECTOR  centeri = this->domains[list[j]].center;
      double x = centeri.x-v.x;
      double y = centeri.y-v.y;
      double z = centeri.z-v.z;
      nearestImage_fast(&x, &y, &z);

      double r2 = x*x + y*y + z*z;
      if (r2 < r2min)
      {
         r2min = r2;
         jmin  = j;
      }
   }

   return list[jmin];
}

void cubic(DDC* ddc)
{
   int lx = ddc->lx;
   int ly = ddc->ly;
   int lz = ddc->lz;
   if (lx*ly*lz != ddc->size && getRank(0) == 0)
   {
      printf("Error in cubic.  Node count mismatch. %d; %d, %d, %d\n",
            ddc->size, lx, ly, lz);
      abortAll(1);
   }

   THREE_VECTOR basis[1] = { {0., 0., 0.} };
   int nBasis = 1;
   generateLattice(ddc, nBasis, basis);
}



void bcc(DDC* ddc)
{
   int lx = ddc->lx;
   int ly = ddc->ly;
   int lz = ddc->lz;
   if (lx*ly*lz*2 != ddc->size && getRank(0) == 0)
   {
      printf("Error in bcc.  Node count mismatch. \n"
            "nTasks=%d; lx, ly, lz=%d, %d, %d\n",
            ddc->size, lx, ly, lz);
      abortAll(1);
   }

   THREE_VECTOR basis[2] = { {-0.25, -0.25, -0.25}, {0.25, 0.25, 0.25} };
   int nBasis = 2;
   generateLattice(ddc, nBasis, basis);
}


void fcc(DDC* ddc)
{
   int lx = ddc->lx;
   int ly = ddc->ly;
   int lz = ddc->lz;
   if (lx*ly*lz*4 != ddc->size && getRank(0) == 0)
   {
      printf("Error in fcc.  Node count mismatch. \n"
            "nTasks=%d; lx, ly, lz=%d, %d, %d\n",
            ddc->size, lx, ly, lz);
      abortAll(1);
   }

   THREE_VECTOR basis[4] = { {-0.25, -0.25, -0.25}, {-0.25,  0.25,  0.25},
      { 0.25, -0.25,  0.25}, { 0.25,  0.25, -0.25} };
   int nBasis = 4;
   generateLattice(ddc, nBasis, basis);
}


void hexagonal(DDC*ddc)
{
   int lx = ddc->lx;
   int ly = ddc->ly;
   int lz = ddc->lz;
   if (lx*ly*lz*2 != ddc->size && getRank(0) == 0)
   {
      printf("HEXAGONAL %d %d %d %d %d\n", ddc->size, 2*lx*ly*lz, lx,ly,lz);
      abortAll(1);
   }

   THREE_VECTOR basis[2] = { {0, -0.25, -0.25}, {0, 0.25, 0.25} };
   int nBasis = 2;
   generateLattice(ddc, nBasis, basis);
}



/** Generates domain centers with an arbitrary basis on a (possibly)
 *  non-uniform cubic lattice.  The lattice may be non-uniform since the
 *  spacing of the lattice planes is determined by the locations of the
 *  centers specified by ddc->[dx, dy, dz].
 *
 *  This routine does not guarantee that all points given by the basis
 *  lie inside their intended cubic cell.  (This is due to the fact that
 *  the domain "center" doesn't necessarily lie at the center of the
 *  cell.  Hence, the width can't be trusted to be equally distributed
 *  on each side of the center.
 *
 *  This routine does not guarantee that all domain centers are located
 *  inside the compuational cell.  (For the same reason basis points
 *  might not be intheir intended cell.)
 *
 *  Note that the basis vectors are added to the center of each cubic
 *  lattice cell.  Hence, for best results you probably want to choose
 *  your lattice so the center of gravity is at (0, 0, 0).
 */
void generateLattice(DDC* ddc, int nBasis, THREE_VECTOR* basis)
{
   assert( ddc->dx != NULL &&  ddc->dy != NULL && ddc->dz != NULL );

   int nTasks = ddc->size;
   int lx = ddc->lx;
   int ly = ddc->ly;
   int lz = ddc->lz;
   if (lx*ly*lz*nBasis != nTasks)
   {
      if (getRank(0) == 0)
         printf("Error in generateLattice.  Node count mismatch.\n"
               "nTasks=%d; lx, ly, lz, nBasis=%d, %d, %d, %d\n",
               nTasks, lx, ly, lz, nBasis);
      abortAll(1);
   }

   double* dx = ddc->dx;
   double* dy = ddc->dy;
   double* dz = ddc->dz;

   double xBinWidth[lx];
   double yBinWidth[ly];
   double zBinWidth[lz];
   int bnd = box_get_boundary_conditions(NULL); 
   mkBinWidth(lx, dx, xBinWidth, bnd & 1);
   mkBinWidth(ly, dy, yBinWidth, bnd & 2);
   mkBinWidth(lz, dz, zBinWidth, bnd & 4);

   THREE_MATRIX h = box_get_h(NULL); 
   THREE_VECTOR corner = box_get_corner(NULL); 
   int m=0;


   for (int ii=0; ii<lx; ++ii)
   {
      THREE_VECTOR ww;
      ww.x = xBinWidth[ii];
      for (int jj=0; jj<ly; ++jj)
      {
         ww.y = yBinWidth[jj];
         for (int kk=0; kk<lz; ++kk)
         {
            ww.z = zBinWidth[kk];
            for (int bb=0; bb<nBasis; ++bb)
            {
               THREE_VECTOR sv;
               sv.x = dx[ii] + basis[bb].x*ww.x;
               sv.y = dy[jj] + basis[bb].y*ww.y;
               sv.z = dz[kk] + basis[bb].z*ww.z;

               THREE_VECTOR center = matrix_vector(h, sv);

               VOP1(center, +=, corner);
               /*
                  ddc->domains.domains[m].radius = 0;
                  ddc->domains.domains[m++].center = center;
                  */
               if(m == (long long int) ddc->domains.local_id) {
                  domainset_setLocalRadius(&ddc->domains,0.0);
                  domainset_setLocalCenter(&ddc->domains,center);
               }
               m++;
            }
         }
      }
   }
   domainset_allGather(&ddc->domains);
}

/** Right now this function does nothing.  Eventually it should take
 *  care of populating random centers.  However, that will have to wait
 *  until domain_set is called only at initialization of the domain
 *  centers.  For now, the random centers are initialized elsewhere as
 *  they have been since Jean-Luc implemented the feature.*/
void nonlattice(DDC* ddc)
{
}

/* This routine takes positions generated as centers of gravity
   of atoms in each processor, after recursive bisection
   subdivision.
   */
void recursive_bisection_domset(DDC *ddc) {
   STATE *state_p = system_getState(NULL);

   assert(state_p != NULL);
   assert(state_p->domcen_active == 1);
   assert(state_p->domcen_np == (int) ddc->domains.size);
   assert(state_p->domcen_pid == (int) ddc->domains.local_id);

   /* Copy data stored in state_p by recursive bisection
      particle distribution into domanset */ {
         /*const*/ /*double (*domcen)[3] = state_p->domcen_vec;
                     const int n = state_p->domcen_np;
                     int i;
                     for(i = 0; i<n; i++) {
                     ddc->domains.domains[i].radius = 0;
                     ddc->domains.domains[i].center.x = domcen[i][0];
                     ddc->domains.domains[i].center.y = domcen[i][1];
                     ddc->domains.domains[i].center.z = domcen[i][2];
                     }
                     */
         {
            const int i = ddc->domains.local_id;
            const double *cp = state_p->domcen_vec[i];
            THREE_VECTOR c;
            c.x = cp[0];
            c.y = cp[1];
            c.z = cp[2];
            domainset_setLocalRadius(&ddc->domains,0.0);
            domainset_setLocalCenter(&ddc->domains,c);
         }
      }

      domainset_allGather(&ddc->domains);
}

double my_rand(void)
{
   return drand48(); 
} 

static int comparePositionsX(const void* av, const void* bv)
{
   const DOMAINX* a = (const DOMAINX*)av;
   const DOMAINX* b = (const DOMAINX*)bv;

   const double aa=a->center.x;
   const double bb=b->center.x;

   if( aa<bb )return -1;
   if( aa>bb )return 1;
   return 0;
}

static int comparePositionsY(const void* av, const void* bv)
{
   const DOMAINX* a = (const DOMAINX*)av;
   const DOMAINX* b = (const DOMAINX*)bv;

   const double aa=a->center.y;
   const double bb=b->center.y;

   if( aa<bb )return -1;
   if( aa>bb )return 1;
   return 0;
}

static int comparePositionsZ(const void* av, const void* bv)
{
   const DOMAINX* a = (const DOMAINX*)av;
   const DOMAINX* b = (const DOMAINX*)bv;

   const double aa=a->center.z;
   const double bb=b->center.z;

   if( aa<bb )return -1;
   if( aa>bb )return 1;
   return 0;
}

void domainset_randomCenters(DOMAINSET* mydomains, int lx, int ly, int lz)
{
   timestamp("Start domainset_randomCenters");

   int nTasks = mydomains->size;


   // generate centers in domain [0,1]**3
   for (int id=0;id<nTasks;id++)
   {
      THREE_VECTOR sv;
      if( lx > 1 )sv.x = (double)my_rand();
      else        sv.x = 0.5;
      if( ly > 1 )sv.y = (double)my_rand();
      else        sv.y = 0.5;
      if( lz > 1 )sv.z = (double)my_rand();
      else        sv.z = 0.5;

      mydomains->domains[id].center.x = sv.x;
      mydomains->domains[id].center.y = sv.y;
      mydomains->domains[id].center.z = sv.z;
   }

   // sort centers, first according to x, then y in bins of sizes ly*lz,
   // then z in bins of sizes lz
   timestamp("Sort random centers");
   if( lx > 1 )
      qsort( mydomains->domains, nTasks, sizeof(DOMAINX), comparePositionsX);
   for(int i=0;i<lx;i++)
      qsort( mydomains->domains+i*ly*lz, ly*lz, sizeof(DOMAINX), comparePositionsY);
   for(int i=0;i<lx*ly;i++)
      qsort( mydomains->domains+i*lz, lz, sizeof(DOMAINX), comparePositionsZ);

   // now put centers in actual computational domain
   THREE_MATRIX h = box_get_h((BOX_STRUCT *)NULL);
   THREE_VECTOR corner = box_get_corner((BOX_STRUCT *)NULL);
   for (int id=0;id<nTasks;id++)
   {
      THREE_VECTOR sv     = mydomains->domains[id].center;
      THREE_VECTOR center = matrix_vector(h, sv);
      VOP1(center, +=, corner);
      mydomains->domains[id].center = center;
   }
   domainset_setLocalCenter(mydomains,mydomains->domains[mydomains->local_id].center); 
}


/**
 *  For a set of nr bin "centers" dr in the range [0, 1] find the
 *  width of each bin.  The bin edges lie half-way between each of the
 *  centers.  Note that the word center is perhaps poorly chosen since
 *  if the centers are not equally spaced the center positions will be
 *  offset from the actual centers of the bins.
 *
 *  Do not assume that half of the width is on each side of the center
 *  of the bin.
 *
 *  This function assumes the dr's are sorted in ascending order. 
 */
static void mkBinWidth(int nr, double* dr, double* binWidth, int periodic)
{
   if (nr == 1)
   {
      binWidth[0] = 1;
      return;
   }

   for (int ii=1; ii<nr-1; ++ii)
      binWidth[ii] = 0.5*(dr[ii+1] - dr[ii-1]);

   if (periodic == 1)
   {
      binWidth[0]    = 0.5*(dr[1] - dr[nr-1] + 1.0);
      binWidth[nr-1] = 0.5*(dr[0] - dr[nr-2] + 1.0);
   }
   else
   {
      binWidth[0]    = 0.5*(dr[0] + dr[1]);
      binWidth[nr-1] = 1 - 0.5*(dr[nr-1] + dr[nr-2]);
   }
}

//unused
/* void domain_print(DOMAINX domain) */
/* { */
/*    printf("DOMAIN: center=(%le,%le,%le), radius=%le\n", */
/*           domain.center.x, domain.center.y, domain.center.z, */
/*           domain.radius); */
/*    fflush(stdout); */
/* } */


/** Returns the domain id of the domain whose center is closest to the
 *  given coordinate */
int domainset_nearestDomain(DOMAINSET* this, THREE_VECTOR r)
{
   const DOMAINX* d = this->domains;
   double r2Min = DIFFSQ(d[0].center, r);
   int r2i = 0;
   for (int ii=1; ii<this->size; ++ii)
   {
      double r2 = DIFFSQ(d[ii].center, r);
      if (r2 < r2Min)
      {
         r2Min = r2;
         r2i = ii;
      }
   }
   return r2i;
}

THREE_VECTOR domainset_getLocalCenter(DOMAINSET* this)
{
   return this->local_domain->center; //this->domains[this->local_id].center;
}

THREE_VECTOR* domainset_getPointerLocalCenter(DOMAINSET* this)
{
   return &(this->local_domain->center);//&this->domains[this->local_id].center;
}

THREE_VECTOR domainset_getRemoteCenter(DOMAINSET* this, int remote_id)
{
   return this->domains[remote_id].center;
}

void domainset_setLocalCenter(DOMAINSET* this, THREE_VECTOR new_center)
{
   this->local_domain->center = new_center;
}

double domainset_getLocalRadius(DOMAINSET* this)
{
   return this->local_domain->radius;//this->domains[this->local_id].radius;
}

void domainset_setLocalRadius(DOMAINSET* this, double radius)
{
   this->local_domain->radius = radius;
}

THREE_VECTOR domainset_getLocalExtent(DOMAINSET* this)
{
   return this->local_domain->extent;//this->domains[this->local_id].radius;
}

void domainset_setLocalExtent(DOMAINSET* this, THREE_VECTOR extent)
{
   this->local_domain->extent = extent;
}
void domainset_setLocalShape(DOMAINSET* this, int shape)
{
   this->local_domain->shape = shape;
}


void domainset_setRadius(DOMAINSET* this, int id, double radius)
{
   this->domains[id].radius = radius;
}

int domainset_overlapBox(DOMAINSET* this, int remote_id, double rcut)
{
   DOMAINX d1 = this->local_domain[0];//domains[this->local_id];
   DOMAINX d2 = this->domains[remote_id];

   THREE_VECTOR v;
   v.x = d1.center.x - d2.center.x;
   v.y = d1.center.y - d2.center.y;
   v.z = d1.center.z - d2.center.z;
   nearestImage(&v.x, &v.y, &v.z);
   double r2 = v.x*v.x + v.y*v.y + v.z*v.z;
   double R = (d1.radius + d2.radius + rcut);
   if (r2 > (R*R) ) return 0; 
   if  ((d1.extent.x + d2.extent.x) + rcut < fabs(v.x)) return 0;
   if  ((d1.extent.y + d2.extent.y) + rcut < fabs(v.y)) return 0;
   if  ((d1.extent.z + d2.extent.z) + rcut < fabs(v.z)) return 0;
   return 1;
}
int domainset_overlapSphere(DOMAINSET* this, int remote_id, double rcut)
{
   DOMAINX d1 = this->local_domain[0];//domains[this->local_id];
   DOMAINX d2 = this->domains[remote_id];

   THREE_VECTOR v;
   v.x = d1.center.x - d2.center.x;
   v.y = d1.center.y - d2.center.y;
   v.z = d1.center.z - d2.center.z;
   nearestImage(&v.x, &v.y, &v.z);
   double r2 = v.x*v.x + v.y*v.y + v.z*v.z;
   double R = (d1.radius + d2.radius + rcut);
   if (r2 > (R*R) ) return 0; 
   return 1;
}
int sphereBoxOverlap(double r, THREE_VECTOR R, THREE_VECTOR Extent)
{
   double *Rptr= &R.x;
   double *Eptr= &Extent.x;
   double d = 0.0;
   for (int i=0;i<3;i++)
   {
      double e = MAX(Rptr[i]-0.5*Eptr[i],0) + MAX(0.5*Eptr[i]-Rptr[i],0);
      if (e <= r) return 0;
      d+=e*e;
   }
   if (d <=r*r) return 1;
   return 0;
}
int domainset_overlap(DOMAINSET* this, int remote_id, double rcut)
{
   enum SHAPE_SHAPE {SphereSphere,SphereBrick,BrickSphere,BrickBrick};
   DOMAINX d1 = this->local_domain[0];//domains[this->local_id];
   DOMAINX d2 = this->domains[remote_id];

   THREE_VECTOR v;
   v.x = d1.center.x - d2.center.x;
   v.y = d1.center.y - d2.center.y;
   v.z = d1.center.z - d2.center.z;
   nearestImage(&v.x, &v.y, &v.z);
   int type = 2*d1.shape + d2.shape;
   //type = SphereSphere;
   switch(type)
   {
      double r2,R;
      case SphereSphere:
      r2 = v.x*v.x + v.y*v.y + v.z*v.z;
      R = (d1.radius + d2.radius) + rcut;
      if (r2 < (R*R)) return 1;
      return 0;
      case SphereBrick:
      return  sphereBoxOverlap(d1.radius+rcut,v,d2.extent);
      case BrickSphere:
      return  sphereBoxOverlap(d2.radius+rcut,v,d1.extent);
      case BrickBrick:
      if ( (d1.extent.x + d2.extent.x + rcut >= fabs(v.x)) &&
            (d1.extent.y + d2.extent.y + rcut >= fabs(v.y)) &&
            (d1.extent.z + d2.extent.z + rcut >= fabs(v.z)) ) return 1;
      return 0;
   }
   return 0;
}

int domain_side(DOMAINX* d1, DOMAINX* d2)
{


   THREE_VECTOR v; 
   v.x = d1->center.x - d2->center.x;
   v.y = d1->center.y - d2->center.y;
   v.z = d1->center.z - d2->center.z;
   backInBox(&v.x, &v.y, &v.z);

   double tol = 1e-5; 
   int x=0;
   int y=0;
   int z=0;
   if (v.x > tol ) x = 1;   if (v.x <-tol ) x =-1;   
   if (v.y > tol ) y = 1;   if (v.y <-tol ) y =-1;   
   if (v.z > tol ) z = 1;   if (v.z <-tol ) z =-1;   
   int l = x + 3*(y + 3*z); 
   if (l > 0) return  1; 
   if (l < 0) return -1; 
   return 0; 
}
void domainset_allGather(DOMAINSET* this)
{
   profile(COMMALLG, START);
   DOMAINX localDomain = *this->local_domain;
   MPI_Allgather(&localDomain,  1, domainx_MPIType(), this->domains, 1, domainx_MPIType(), COMM_LOCAL);
   profile(COMMALLG, END);
}


double domainset_getMaxRadius(DOMAINSET* this)
{
   double max_radius =0.0;
   for (int i=0;i<this->size;i++)
      max_radius = MAX( max_radius,this->domains[i].radius) ; 
   return max_radius;
}

double domainset_setRadiusToMax(DOMAINSET* this)
{
   double max_radius =0.0;
   MPI_Allreduce(&(this->domains[this->local_id].radius), &max_radius, 1, MPI_DOUBLE, MPI_MAX, COMM_LOCAL);
   for (int i=0;i<this->size;i++)
      this->domains[i].radius = max_radius;
   return max_radius;
}


/** In May 2011 I removed all the code from cubic, fcc, etc, that found
 *  the bounding box for the local domain.  This greatly simplifies
 *  domain_set (especially in non-lattice cases), but breaks the
 *  routingHack code in pppm_makeRho.c.  However, I did not want to throw
 *  away the rouingHack code just yet.  It will still be useful if there
 *  was a way to get a handle on the bounding boxes.  (It also works
 *  much better in large systems.)  So, this is a placeholder for a
 *  function that is able to look at the current domain centers and find
 *  the boundingbox of the local domain.  If you're thinking of writing
 *  one you should probably plan to handle the case where the centers
 *  aren't on a lattice because the load balancing code has been moving
 *  them around.  Older versions of computeLocalVoronoiRadius may be
 *  helpful as a start. - dfr
 */
void domain_getBoundingBox(THREE_VECTOR* rMin, THREE_VECTOR* rMax)
{
   if (getRank(0) != 0) return;
   printf("Game Over\n"
         "Don't try to use routingHack=1 in PPPM.\n"
         "It doesn't work anymore.\n"
         "Blame Dave.\n"
         "If you really need it, maybe he will fix it for you.\n");
   assert(1==0);
}


/* Local Variables: */
/* tab-width: 3 */

/* End: */
