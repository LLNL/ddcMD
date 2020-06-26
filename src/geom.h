#ifndef GEOM_H
#define GEOM_H

#include <stdlib.h>
#include "three_algebra.h"
#include "box.h"
#include "gid.h"

enum GEOM_METHOD {GEOM_NONE, GEOM_DEFAULT, GEOM_FAST, GEOM_P7N1ORTHO };

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))

typedef struct geombox_st
{
	THREE_VECTOR r;
	unsigned ix,iy,iz;
	unsigned int key, index, type;
	int nn, first, firstlocal;
	int *nlist;
	THREE_VECTOR **R;
} GEOMBOX;

/** A PARTICLESET is a way to abstract the representation of a
 *  collection of particle data so that the coordinates, gids, and types
 *  can be accessed without knowing any of the details of the classes in
 *  which the data is actually stored.  This allows GEOM to directly
 *  access particle data in nearly any form and allows us to avoid
 *  making a copy in a GEOM specific representation.
 *
 *  No actual data is stored here, only pointers to the data.
 *
 *  To use a PARTICLESET, the positions, x, y, and z, must be stored as
 *  doubles, the gids must be gid_type, and the types must be ints.  The
 *  data storage must also be arranged such that each field of the data
 *  for two successive particles are separated by a fixed number of
 *  bytes.  For example, if the the x-coordinates are stored in an array
 *  of doubles, the stride is sizeof(double).  If the x-, y-, and
 *  z-coordinates are stored as an array of THREE_VECTORs the stride is
 *  sizeof(THREE_VECTOR).  In either case
 *               *(double*)((char*)xtpr+ii*stride)
 *  is the value of the xcoordinate of the iith particle.
 *
 *  We expect that the members of the stuct will be initialized such that:
 *
 *  - xptr, yptr, zptr, global_index, and type all point to their
 *    respective data arrays for particle 0.
 *  - stridedouble is the number of bytes between to successive
 *    values of x, y, or z.
 *  - strideint is the number of of bytes between to successive
 *    values of type.
 *  - stridelong64 is the number of of bytes between to successive
 *    values of the gid.
 *  - number_particles and number_local are pointers to the global and
 *    local number of particles (respectively).
 *  - center points to a THREE_VECTOR that specifies a
 *    coordinate that is somewhere near the center of mass for the
 *    particles.  This coordinate will be used as the origin of the
 *    coordinate system for the particle on the GEOM grid.  This ensure
 *    that the particles have coordinates within a compact region of
 *    space. 
 *  
 */
typedef struct set_str
{
   unsigned *number_particles, *number_local;
   int stridedouble, strideint, stridelong64;
	double *xptr, *yptr, *zptr;
	gid_type *global_index;
	int *type;
	const THREE_VECTOR* center;
} PARTICLESET;

typedef struct geom_pinfo_st
{
	GEOMBOX* box;  // which cell is this particle in
	int next;      // next particle in same cell as this one.
} GEOM_PINFO;


/**
 *  GEOM is used to find pairs of particles.  It uses the standard trick
 *  of gridding space and only searching a subset of grid cells to
 *  accelerate the pair finding.
 *
 *  Some members of this struct are unused.  They are dinosaur tracks,
 *  or remnants of code that was never completed.  These include
 *  boundary_condition, small, and box_nlist_update.
 */
typedef struct geom_str
{
	enum GEOM_METHOD method;
   BOX_STRUCT *compBox; 
	THREE_MATRIX *h, *hinv; // pointers to data in BOX
   THREE_VECTOR rCorner; 
	THREE_VECTOR d;
	THREE_VECTOR min,max; 
	THREE_VECTOR cellSpan;
	THREE_VECTOR n0, n1, n2;
	int *boundary_condition; // pointer to data in BOX
	double minBoxSide;
	double rcut; 
	int small, nx, ny, nz, nbox;
	int box_nlist_update; 
	int boxCapacity;
	GEOMBOX *box;
	int pinfoCapacity;
	GEOM_PINFO* pinfo;
	THREE_VECTOR* r;
	int boxOffsetCapacity;
	int* boxOffset;
} GEOM;

GEOM* geom_init(double rCut, BOX_STRUCT *box);
void  geom_free(GEOM* self);
void GeomBox(GEOM* geom, PARTICLESET* ps);
void geomMethodSet(GEOM* self, enum GEOM_METHOD method);
PARTICLESET* ParticleSet(
	PARTICLESET* p,
	double *xptr, double *yptr, double *zptr,	int stridedouble,
   gid_type* global_index, int stridelong64,
	int* type, int strideint,
   unsigned* number_particles, unsigned* number_local,
   const THREE_VECTOR* center);

#endif

/* Local Variables: */
/* tab-width: 3 */
/* End: */
