#include "preduce.h"

#include <math.h>
#include <stdlib.h> // NULL
#include <assert.h>

#define XX   0
#define YX   1
#define ZX   2
#define XY   3
#define YY   4
#define ZY   5
#define XZ   6
#define YZ   7
#define ZZ   8
#ifdef ANSI_F77
#define RINT rint
#endif
#ifdef AIX
#define RINT nearest
#endif
#ifdef HPUX
#define RINT nint
#endif
#ifndef RINT
#define RINT rint
#endif
double nearbyint(double);
static double da, db, dc, dx, dy, dz;
static double* h = NULL;
static double* hi = NULL;
static double hhxx,hhyy,hhzz; 
static double hxx,hyy,hzz; 
static int bndcdn=-1;


static int* _nShifts = NULL;
static THREE_INT** _lShift= NULL;
static double* _r2Safe=NULL;


void Preduce(double *x, double *y, double *z);
void Preduce_image(double *x, double *y, double *z);

static int orthorhombicBox(THREE_MATRIX* h);

void (*nearestImage)(double*, double*, double*) = Preduce_image;
void (*backInBox)(double*, double*, double*) = Preduce;
void (*nearestImage_fast)(double*, double*, double*) = Preduce_image;
void (*backInBox_fast)(double*, double*, double*) = Preduce;


void Pset(double *hptr, double *hiptr, int *bndptr, double* r2Safe, int* nShifts, THREE_INT** lShift)
{
	h = hptr;
	hi = hiptr;
	bndcdn = *bndptr;
	hhxx = 0.5*h[XX]; 
	hhyy = 0.5*h[YY]; 
	hhzz = 0.5*h[ZZ]; 
	hxx = h[XX]; 
	hyy = h[YY]; 
	hzz = h[ZZ];
	_r2Safe = r2Safe;
	_nShifts = nShifts;
	_lShift = lShift;

	// set preduce methods
	PsetMethod(*bndptr, (THREE_MATRIX*) hptr);
}

	
void dpreduce(double *x, double *y, double *z)
{
	switch (bndcdn)
	{
	  case 0:
	   break;
	case 1:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		*x = h[XX]*da;
		*y = h[XY]*da;
		*z = h[XZ]*da;
		break;
	case 2:
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		*x = h[YX]*db;
		*y = h[YY]*db;
		*z = h[YZ]*db;
		break;
	case 3:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		*x = h[XX]*da + h[YX]*db;
		*y = h[XY]*da + h[YY]*db;
		*z = h[XZ]*da + h[YZ]*db;
		break;
	case 4:
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x = h[ZX]*dc;
		*y = h[ZY]*dc;
		*z = h[ZZ]*dc;
		break;
	case 5:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x = h[XX]*da + h[ZX]*dc;
		*y = h[XY]*da + h[ZY]*dc;
		*z = h[XZ]*da + h[ZZ]*dc;
		break;
	case 6:
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x = h[YX]*db + h[ZX]*dc;
		*y = h[YY]*db + h[ZY]*dc;
		*z = h[YZ]*db + h[ZZ]*dc;
		break;
	case 7:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x = h[XX]*da + h[YX]*db + h[ZX]*dc;
		*y = h[XY]*da + h[YY]*db + h[ZY]*dc;
		*z = h[XZ]*da + h[YZ]*db + h[ZZ]*dc;
		break;
	  default:
	   assert(1==0);
	}
}
void PreduceB0(double *x, double *y, double *z)
{
	return; 
}

void PreduceOrthorhombicB7_OneLatticeReductionNEW(double *x, double *y, double *z)
{

		if (*x > hhxx ) {  *x -= hxx; }
		if (*x < -hhxx ) { *x += hxx; }

		if (*y > hhyy )  { *y -= hyy; }
		if (*y < -hhyy ) { *y += hyy; }

		if (*z > hhzz ) {  *z -= hzz; }
		if (*z < -hhzz ) { *z += hzz; }
}
void PreduceOrthorhombicB7_OneLatticeReduction(double *x, double *y, double *z)
{
		double hxx,hyy,hzz; 
		hxx = h[XX]; hyy = h[YY]; hzz = h[ZZ]; 

		if (*x > 0.5*hxx ) {  *x += dx=-hxx; }
		if (*x < -0.5*hxx ) { *x += dx= hxx; }

		if (*y > 0.5*hyy )  { *y += dy=-hyy; }
		if (*y < -0.5*hyy ) { *y += dy= hyy; }

		if (*z > 0.5*hzz ) {  *z += dz=-hzz; }
		if (*z < -0.5*hzz ) { *z += dz= hzz; }
}
void VPreduceOrthorhombicB7_OneLatticeReduction(double *x, double *y, double *z, unsigned  n)
{
		double hxx,hyy,hzz; 
		hxx = h[XX]; hyy = h[YY]; hzz = h[ZZ]; 
		double hhxx=0.5*hxx; 
		double hhyy=0.5*hyy; 
		double hhzz=0.5*hzz; 

		for (unsigned i=0;i<n;i++)
		{
		if (x[i] > hhxx ) { x[i] -=  hxx; }
		if (x[i] <-hhxx ) { x[i] +=  hxx; }

		if (y[i] > hhyy ) { y[i] -=  hyy; }
		if (y[i] <-hhyy ) { y[i] +=  hyy; }

		if (z[i] > hhzz ) { z[i] -=  hzz; }
		if (z[i] <-hhzz ) { z[i] +=  hzz; }
		}
}
void PreduceOrthorhombicB6_OneLatticeReduction(double *x, double *y, double *z)
{
		double hyy,hzz; 
		hyy = h[YY]; hzz = h[ZZ]; 

		if (*y > 0.5*hyy ) { *y += dy = -hyy ; }
		if (*y < -0.5*hyy ) { *y += dy = hyy; }

		if (*z > 0.5*hzz ) {  *z += dz =-hzz ; }
		if (*z < -0.5*hzz ) { *z += dz = hzz; }
}
void PreduceOrthorhombicB5_OneLatticeReduction(double *x, double *y, double *z)
{
		double hxx,hzz; 
		hxx = h[XX]; hzz = h[ZZ]; 

		if (*x > 0.5*hxx ) {  *x += dx =-hxx ; }
		if (*x < -0.5*hxx ) { *x += dx = hxx; }

		if (*z > 0.5*hzz ) {  *z += dz =-hzz ; }
		if (*z < -0.5*hzz ) { *z += dz = hzz; }
}
void PreduceOrthorhombicB4_OneLatticeReduction(double *x, double *y, double *z)
{
		double hzz; 
		hzz = h[ZZ]; 

		if (*z > 0.5*hzz ) {  *z += dz =-hzz ; }
		if (*z < -0.5*hzz ) { *z += dz = hzz; }
}
void PreduceOrthorhombicB3_OneLatticeReduction(double *x, double *y, double *z)
{
		double hxx,hyy; 
		hxx = h[XX]; hyy = h[YY]; 

		if (*x > 0.5*hxx ) {  *x += dx =-hxx ; }
		if (*x < -0.5*hxx ) { *x += dx = hxx; }

		if (*y > 0.5*hyy ) { *y += dy = -hyy ; }
		if (*y < -0.5*hyy ) { *y += dy = hyy; }

}
void PreduceOrthorhombicB2_OneLatticeReduction(double *x, double *y, double *z)
{
		double hyy; 
		hyy = h[YY]; 

		if (*y > 0.5*hyy ) { *y += dy = -hyy ; }
		if (*y < -0.5*hyy ) { *y += dy = hyy; }

}
void PreduceOrthorhombicB1_OneLatticeReduction(double *x, double *y, double *z)
{
		double hxx = h[XX]; 
		if (*x > 0.5*hxx ) *x += -hxx ;
		if (*x < -0.5*hxx ) *x += hxx;
}


/** This function finds a nearest the nearest image of the vector (x, y,
 * z).  If the vector found by a back-in-box operation is within the
 * inscribed sphere, then we can't possibly do any better.  Otherwise,
 * we search over all lattice shifts that possibly intersect the
 * Wigner-Seitz cell. */
void Preduce_image(double *x, double *y, double *z)
{
   Preduce(x, y, z);

   if ((*x)*(*x) + (*y)*(*y) + (*z)*(*z) < *_r2Safe)
      return;
   
   double xMin = *x;
   double yMin = *y;
   double zMin = *z;
   double r2Min = xMin*xMin + yMin*yMin + zMin*zMin;
   for (int ii=0; ii<*_nShifts; ++ii)
   {
      int ix = (*_lShift)[ii].x;
      int iy = (*_lShift)[ii].y;
      int iz = (*_lShift)[ii].z;
      
      double xTmp = *x + h[XX]*ix + h[YX]*iy + h[ZX]*iz;
      double yTmp = *y + h[XY]*ix + h[YY]*iy + h[ZY]*iz;
      double zTmp = *z + h[XZ]*ix + h[YZ]*iy + h[ZZ]*iz;
      double r2Tmp = xTmp*xTmp + yTmp*yTmp + zTmp*zTmp;
      if (r2Tmp < r2Min)
      {
	 r2Min = r2Tmp;
	 xMin = xTmp;
	 yMin = yTmp;
	 zMin = zTmp;
      }
   }
   *x = xMin;
   *y = yMin;
   *z = zMin;
}




void Preduce(double *x, double *y, double *z)
{
	switch (bndcdn)
	{
	  case 0:
	   break;
	case 1:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		*x += dx = h[XX]*da;
		*y += dy = h[XY]*da;
		*z += dz = h[XZ]*da;
		break;
	case 2:
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		*x += dx = h[YX]*db;
		*y += dy = h[YY]*db;
		*z += dz = h[YZ]*db;
		break;
	case 3:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		*x += dx = h[XX]*da + h[YX]*db;
		*y += dy = h[XY]*da + h[YY]*db;
		*z += dz = h[XZ]*da + h[YZ]*db;
		break;
	case 4:
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x += dx = h[ZX]*dc;
		*y += dy = h[ZY]*dc;
		*z += dz = h[ZZ]*dc;
		break;
	case 5:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x += dx = h[XX]*da + h[ZX]*dc;
		*y += dy = h[XY]*da + h[ZY]*dc;
		*z += dz = h[XZ]*da + h[ZZ]*dc;
		break;
	case 6:
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x += dx = h[YX]*db + h[ZX]*dc;
		*y += dy = h[YY]*db + h[ZY]*dc;
		*z += dz = h[YZ]*db + h[ZZ]*dc;
		break;
	case 7:
		da = -RINT(hi[XX]*(*x) + hi[YX]*(*y) + hi[ZX]*(*z));
		db = -RINT(hi[XY]*(*x) + hi[YY]*(*y) + hi[ZY]*(*z));
		dc = -RINT(hi[XZ]*(*x) + hi[YZ]*(*y) + hi[ZZ]*(*z));
		*x += dx = h[XX]*da + h[YX]*db + h[ZX]*dc;
		*y += dy = h[XY]*da + h[YY]*db + h[ZY]*dc;
		*z += dz = h[XZ]*da + h[YZ]*db + h[ZZ]*dc;
		break;
	  default:
	   assert(1==0);
	}
}

void preduceGetR(THREE_VECTOR*R)
{
	R->x = dx;
	R->y = dy;
	R->z = dz;
}
void preduceGetS(THREE_VECTOR*S)
{
	
	S->x = da;
	S->y = db;
	S->z = dc;
}

void Preduce1(double *x, double *y, double *z, double *daptr, double *dbptr, double *dcptr)
{
	Preduce(x, y, z);
	*daptr = da;
	*dbptr = db;
	*dcptr = dc;
}

void Preduce2(double *x, double *y, double *z, double *dxptr, double *dyptr, double *dzptr)
{
	Preduce(x, y, z);
	*dxptr = dx;
	*dyptr = dy;
	*dzptr = dz;
}

void preduce1(double *x, double *y, double *z, double *daptr, double *dbptr, double *dcptr)
{
	Preduce1(x, y, z, daptr, dbptr, dcptr);
}

void PreduceUnitBox(THREE_VECTOR*v)
{
	switch (bndcdn)
	{
	  case 0:
	   break;
	case 1:
		v->x -= floor(v->x);
		break;
	case 2:
		v->y -= floor(v->y);
		break;
	case 3:
		v->x -= floor(v->x);
		v->y -= floor(v->y);
		break;
	case 4:
		v->z -= floor(v->z);
		break;
	case 5:
		v->x -= floor(v->x);
		v->z -= floor(v->z);
		break;
	case 6:
		v->y -= floor(v->y);
		v->z -= floor(v->z);
		break;
	case 7:
		v->x -= floor(v->x);
		v->y -= floor(v->y);
		v->z -= floor(v->z);
		break;
	  default:
	   assert(1==0);
	}
}

void PreduceMinBox(THREE_VECTOR*v)
{
	switch (bndcdn)
	{
	  case 0:
	   break;
	case 1:
		v->x -= rint(v->x);
		break;
	case 2:
		v->y -= rint(v->y);
		break;
	case 3:
		v->x -= rint(v->x);
		v->y -= rint(v->y);
		break;
	case 4:
		v->z -= rint(v->z);
		break;
	case 5:
		v->x -= rint(v->x);
		v->z -= rint(v->z);
		break;
	case 6:
		v->y -= rint(v->y);
		v->z -= rint(v->z);
		break;
	case 7:
		v->x -= rint(v->x);
		v->y -= rint(v->y);
		v->z -= rint(v->z);
		break;
	  default:
	   assert(1==0);
	}
}

void PsetMethod(unsigned pbc, THREE_MATRIX* h)
{
   // most conservative settings
   nearestImage      = Preduce_image;
   nearestImage_fast = Preduce_image;
   backInBox         = Preduce;
   backInBox_fast    = Preduce;
   
   if (pbc==0)
   {
      nearestImage      = PreduceB0;
      nearestImage_fast = PreduceB0;
      backInBox         = PreduceB0; 
      backInBox_fast    = PreduceB0;
   }
   
   if (orthorhombicBox(h) == 1)
   {
      nearestImage = Preduce;
      backInBox    = Preduce;
      void (*method)(double*, double*,  double*) = Preduce;
      if (pbc==7) method = PreduceOrthorhombicB7_OneLatticeReduction; 
      if (pbc==6) method = PreduceOrthorhombicB6_OneLatticeReduction; 
      if (pbc==5) method = PreduceOrthorhombicB5_OneLatticeReduction; 
      if (pbc==4) method = PreduceOrthorhombicB4_OneLatticeReduction; 
      if (pbc==3) method = PreduceOrthorhombicB3_OneLatticeReduction; 
      if (pbc==2) method = PreduceOrthorhombicB2_OneLatticeReduction; 
      if (pbc==1) method = PreduceOrthorhombicB1_OneLatticeReduction;
      backInBox_fast = method;
      nearestImage_fast = method;
   }
}

/** Returns 1 if h is an orthorhombic box, 0 otherwise. */
int orthorhombicBox(THREE_MATRIX* h)
{
   double eps = 1e-10;
   if ( (fabs(h->xy) > eps) ||
	(fabs(h->xz) > eps) ||
	(fabs(h->yx) > eps) ||
	(fabs(h->yz) > eps) ||
	(fabs(h->zx) > eps) ||
	(fabs(h->zy) > eps) )
      return 0;
   return 1;
}




/* Local Variables: */
/* tab-width: 3 */
/* End: */
