#include "box.h"
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "mpiUtils.h"
#include "three_algebra.h"
#include "object.h"
#include "units.h"
#include "eq.h"
#include "mpiUtils.h"
#include "preduce.h"
#include "ddcMalloc.h"
#include "gid.h"

// We provide a prototype so that we don't have to include system.h.
// That allows us to reuse this code for other projects.
struct system_st;
gid_type system_getNglobal(struct system_st* system);


#define MIN(A,B) ((A) < (B) ? (A) : (B))

static BOX_STRUCT *_box_current;
static void misspelledKeywordError(void);
static void updatePreduce(BOX_STRUCT* box);
static void findCellProperties(BOX_STRUCT* box);


// pbc=1 (001): periodic in x
// pbc=2 (010): periodic in y
// pbc=3 (011): periodic in x and y
// pbc=4 (100): periodic in z
// pbc=5 (101): periodic in x and z
// pbc=6 (110): periodic in y and z
// pbc=7 (111): periodic in x, y and z
void updateH0(BOX_STRUCT *box,THREE_MATRIX h)
{
   box->h0=h; 
   box->hfac =  matrix_matrix(h,box->hinv); 
   if (matrix_equal_tol(box->hfac, I_3x3, 1e-14)) box->hfac=I_3x3;
   box->volume = matinv(&box->h0, &box->hinv);
   box->corner = matrix_vector(box->h0, box->reducedcorner);
   updatePreduce(box);
}

BOX_STRUCT *box_init(void *parent, char *name)
{
   BOX_STRUCT *box;
   box = (BOX_STRUCT *) object_initialize(name, "BOX", sizeof(BOX_STRUCT));
   box->parent = parent; 
   OBJECT* obj = (OBJECT*) box;
   object_get(obj, "type", &box->type, STRING, 1, "GENERAL");
   // for backwards compatibility, use old keyword if present.
   if (object_testforkeyword(obj, "bndcdn"))
      object_get(obj, "bndcdn", &box->pbc, INT, 1, "7");
   else
      object_get(obj, "pbc", &box->pbc, INT, 1, "7");
   if (object_testforkeyword(obj, "bdncdn") || object_testforkeyword(obj, "bndcnd"))
      misspelledKeywordError();

   object_get(obj, "h", &(box->h0), WITH_UNITS, 9, "1 0 0 0 1 0 0 0 1","l",NULL);
   object_get(obj, "dhdt", &(box->dhdt), WITH_UNITS, 9, "0 0 0 0 0 0 0 0 0","l",NULL);
   object_get(obj, "reducedcorner", &(box->reducedcorner), DOUBLE, 3, "-0.5 -0.5 -0.5");
   box->volume = matinv(&box->h0, &box->hinv);
   box->corner = matrix_vector(box->h0, box->reducedcorner);
   box->itype=GENERAL; 
   if (strcmp(box->type ,"ORTHORHOMBIC")==0) box->itype=ORTHORHOMBIC; 
   box->time_dependence = BOX_NONE;
   box->bodyDiag = vzero;
   box->r2Inner = 0;
   box->r2Outer = 0;
   box->lvSize = 0;
   box->lvCapacity = 27;
   box->lv = ddcMalloc(box->lvCapacity * sizeof(THREE_INT));
   findCellProperties(box);
   _box_current = box;

   Pset((double*) &box->h0, (double*) &box->hinv, &box->pbc,
         &box->r2Inner, &box->lvSize, &box->lv);

   box->nAffineUpdates = 0;
   box->affineUpdates = NULL;

   return box;
}

void box_write(BOX_STRUCT*box, FILE*file)
{
   char *fmt = "%21.14e %21.14e %21.14e\n      %21.14e %21.14e %21.14e\n      %21.14e %21.14e %21.14e;\n";
   THREE_MATRIX h0 = box->h0;
   THREE_MATRIX dhdt = box->dhdt;
   double convert_factor = units_convert(1.0,NULL,"l"); 
   fprintf(file, "%s %s {\n h  = ", box->name, box->objclass);
   MATSCALE(h0, convert_factor);
   MFPRINT(file, fmt, h0);
   if (!matrix_equal_tol(dhdt, mzero, 1e-15))
   {
   convert_factor = units_convert(1.0,NULL,"l/t"); 
   MATSCALE(dhdt, convert_factor);
   fprintf(file,"dhdt = "); 
   MFPRINT(file, fmt, dhdt);
   }
   fprintf(file, "}\n");
}

void box_put(BOX_STRUCT*box, int cmd, void *ptr)
{
   double volume_new;
   THREE_MATRIX h;
   if (box == NULL) box = _box_current;
   switch (cmd)
   {
      case BOX_TIMESTART:
         box->time = *(double*)ptr; 
         break; 
      case BOX_TIME:
         h = box->h0;
         if (box->time_dependence != BOX_NONE) 
         {
            double newTime = *(double*)ptr; 
            h=boxPrescriptiveTime(box,newTime); 
            box->time = newTime;
         }
         updateH0(box,h); 
         break; 
      case HO:
         h = *(THREE_MATRIX*)ptr; 
         updateH0(box,h); 
         break ; 
      case DHDT:
         box->dhdt = *(THREE_MATRIX*)ptr; 
         break ; 
      case HFAC:
         box->hfac= *(THREE_MATRIX *)ptr; 
         break; 
      case VOLUME:
         h = box->h0; 
         volume_new = *(double *)ptr;
         double scale = cbrt(volume_new/box->volume);
         MATSCALE(h, scale);
         updateH0(box,h); 
         break;
      case RESET:
         box->corner = matrix_vector(box->h0, box->reducedcorner);
         box->volume = matinv(&box->h0, &box->hinv);
         updatePreduce(box);
         break; 
   }
}
BOX_STRUCT * box_getBox(BOX_STRUCT *box)
{
   if ( box == NULL ) box = _box_current ;
   return box; 
}
int box_getType(BOX_STRUCT *box)
{
   if ( box == NULL ) box = _box_current ;
   return box->itype; 
}
THREE_VECTOR box_get_diagonal(BOX_STRUCT *box)
{ 
   if ( box == NULL ) box = _box_current ;
   THREE_VECTOR ur3={1.,1.,1.};
   THREE_VECTOR d = matrix_vector(box->h0, ur3);
   return d;
}

THREE_MATRIX box_get_h(BOX_STRUCT *box)
{
   if ( box == NULL ) box = _box_current ;
   return box->h0;
}
THREE_MATRIX box_get_dhdt(BOX_STRUCT *box)
{
   if ( box == NULL ) box = _box_current ;
   return box->dhdt;
}
THREE_VECTOR box_get_corner(BOX_STRUCT *box) // get lower left corner
{ 
   if ( box == NULL ) box = _box_current;
   return box->corner;
}
THREE_VECTOR box_get_reducedcorner(BOX_STRUCT *box)
{ 
   if ( box == NULL ) box = _box_current;
   return box->reducedcorner;
}
THREE_VECTOR box_get_urcorner(BOX_STRUCT *box)
{ 
   if ( box == NULL ) box = _box_current;
   THREE_VECTOR ur3={1.,1.,1.};
   THREE_VECTOR ur = matrix_vector(box->h0, ur3);
   ur.x+=box->corner.x;
   ur.y+=box->corner.y;
   ur.z+=box->corner.z;

   return ur;
}
double box_get_volume(BOX_STRUCT *box)
{
   if ( box == NULL ) box = _box_current ;
   return box->volume = matinv(&box->h0, &box->hinv);
}
int  box_get_boundary_conditions(BOX_STRUCT *box)
{
   if (box == NULL) box=_box_current;
   return box->pbc;
}
double box_get_minspan(BOX_STRUCT  *box)
{
   THREE_VECTOR lv[] = { {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}, {0, 1, -1}, {1, 0, -1}, {1, -1, 0}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 0} };
   int nlv=13; 
   if ( box == NULL ) box = _box_current ; 
   THREE_VECTOR r = matrix_vector(box->h0, lv[0]);
   double r2min = r.x*r.x + r.y*r.y + r.z*r.z;
   for (int i = 1; i < nlv; i++)
   {
      r = matrix_vector(box->h0, lv[i]);
      double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
      if (r2 < r2min) r2min = r2;
   }
   return sqrt(r2min); 
}
THREE_VECTOR box_get_boundingbox(BOX_STRUCT *box) 
{ 
   THREE_VECTOR r,rmin,rmax, lv[] = { {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1} };
   int i, nlv = 7;
   if ( box == NULL ) box = _box_current ; 
   rmin = rmax = matrix_vector(box->h0, lv[0]);
   for (i = 1; i < nlv;i++)
   {
      r = matrix_vector(box->h0, lv[i]);
      if (r.x > rmax.x ) rmax.x = r.x ;
      if (r.y > rmax.y ) rmax.y = r.y ;
      if (r.z > rmax.z ) rmin.z = r.z ;
      if (r.x < rmin.x ) rmin.x = r.x ;
      if (r.y < rmin.y ) rmin.y = r.y ;
      if (r.z < rmin.z ) rmin.z = r.z ;
   }
   r.x = rmax.x - rmin.x ; 
   r.y = rmax.y - rmin.y ; 
   r.z = rmax.z - rmin.z ; 
   return r; 
}

void box_get(BOX_STRUCT*box, int cmd, void *ptr)
{
   double r2, r2min;
   THREE_VECTOR r, lv[] = { {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}, {0, 1, -1}, {1, 0, -1}, {1, -1, 0}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 0} };
   static THREE_VECTOR R[27];
   int i, nlv = 13;

   if (box == NULL) box = _box_current;
   switch (cmd)
   {
      case VOLUME:
         box->volume = matinv(&box->h0, &box->hinv);
         *(double *)ptr = box->volume;
         break;
      case HO_PTR:
         *(THREE_MATRIX **) ptr = &box->h0;
         break;
      case HINV_PTR:
         *(THREE_MATRIX **) ptr = &box->hinv;
         break;
      case CORNER_PTR:
         *(THREE_VECTOR **) ptr = &box->corner;
         break;
      case PBC_PTR:
         *(int **)ptr = &box->pbc;
         break;
      case MINSPAN:
         r = matrix_vector(box->h0, lv[0]);
         r2min = r.x*r.x + r.y*r.y + r.z*r.z;
         for (i = 1; i < nlv; i++)
         {
            r = matrix_vector(box->h0, lv[i]);
            r2 = r.x*r.x + r.y*r.y + r.z*r.z;
            if (r2 < r2min) r2min = r2;
         }
         *(double *)ptr = sqrt(r2min);
         break;
      case LATTICEVECTORS:

         {
            THREE_VECTOR *r;
            r = (THREE_VECTOR *)ptr; 
            r[0].x = box->h0.xx; 
            r[0].y = box->h0.yx; 
            r[0].z = box->h0.zx; 
            r[1].x = box->h0.xy; 
            r[1].y = box->h0.yy; 
            r[1].z = box->h0.zy; 
            r[2].x = box->h0.xz; 
            r[2].y = box->h0.yz; 
            r[2].z = box->h0.zz; 
         }
         break;
      case RECIP_LATTICEVECTORS:

         {
            THREE_VECTOR *recip;
            recip = (THREE_VECTOR *)ptr; 
            recip[0].x = 2.0*M_PI*box->hinv.xx; 
            recip[0].y = 2.0*M_PI*box->hinv.xy; 
            recip[0].z = 2.0*M_PI*box->hinv.xz; 
            recip[1].x = 2.0*M_PI*box->hinv.yx; 
            recip[1].y = 2.0*M_PI*box->hinv.yy; 
            recip[1].z = 2.0*M_PI*box->hinv.yz; 
            recip[2].x = 2.0*M_PI*box->hinv.zx; 
            recip[2].y = 2.0*M_PI*box->hinv.zy; 
            recip[2].z = 2.0*M_PI*box->hinv.zz; 
         }
         break;
      case NEAREST_IMAGES:
         {
            int ix,iy,iz; 
            THREE_VECTOR s; 
            i=0; 
            for (iz = -1 ; iz  <= 1;iz ++) 
               for (iy = -1 ; iy  <= 1;iy ++) 
                  for (ix = -1 ; ix  <= 1;ix ++) 
                  {
                     s.x = ix;  s.y = iy ; s.z = iz; 
                     R[i++] = matrix_vector(box->h0, s);
                  }
            *(THREE_VECTOR **)ptr = R; 
         }
         break; 
      case HFAC :
         *(THREE_MATRIX*)ptr = box->hfac; 
         break ; 
   }
}

/**
 *  This function updates various properties of the simulation box that
 *  are needed for nearest image calculations in non-orthorhombix boxes.
 *  The idea is that this function should be called each time the box
 *  size or shape changes.  Specifically, we find:
 *
 *  - r2Inner:  squared radius of the largest inscribed sphere.
 *  - r2Outer:  squared radius of the circumsphere.
 *  - bodyDiag: longest body diagonal vector
 *  - lv:       set of integer lattice shift vectors to all image
 *              cells that overlap the Wigner-Seitz cell.  It is
 *              possible that some cells in this list do not overlap
 *              the W_S cell, but all that do touch are in the list.
 *  - lvSize:   number of entries in lv
 *
 *  The current algorithm to find possibly overlapping cells is the
 *  simplest and most robust that I could think of.  We find the
 *  circumsphere of the simulation cell.  The Wigner-Seitz cell must lie
 *  within this sphere.  We now find every image cell whose circumsphere
 *  intersects the circumsphere of the simulation cell.  These image
 *  cells might overlap the W-S cell.  For cases involving high aspect
 *  ratio simulation cells this can be a disasterously large over
 *  estimate of possibly overlapping images.  However we can fix that
 *  problem later.
 *
 *  If you're calling this function you almost certainly want to call
 *  PsetMethod too to update the preduce function pointers (or you could
 *  just call updatePreduce for one stop shopping).
 */
void findCellProperties(BOX_STRUCT* box)
{
   // find longest body-diagonal
   THREE_VECTOR corners[4] = { {1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}, {1.0, -1.0, 1.0}, {1.0, 1.0,-1.0} };

   double maxLen = -1;
   for (unsigned ii=0; ii<4; ++ii)
   {
      THREE_VECTOR ri = matrix_vector(box->h0, corners[ii]);
      double len = dot1(ri, ri);
      if (len > maxLen)
      {
         box->bodyDiag = ri;
         maxLen = len;
      }
   }

   box->r2Outer = maxLen/4.0;

   // compute rInner
   THREE_VECTOR a0 = {1, 0, 0};
   THREE_VECTOR a1 = {0, 1, 0};
   THREE_VECTOR a2 = {0, 0, 1};
   a0 = matrix_vector(box->h0, a0);
   a1 = matrix_vector(box->h0, a1);
   a2 = matrix_vector(box->h0, a2);

   double a0_len = sqrt(dot1(a0, a0));
   double a1_len = sqrt(dot1(a1, a1));
   double a2_len = sqrt(dot1(a2, a2));

   THREE_VECTOR n0 = cross(&a1, &a2);
   THREE_VECTOR n1 = cross(&a2, &a0);
   THREE_VECTOR n2 = cross(&a0, &a1);
   VSCALE(n0, 1.0/(a1_len * a2_len));
   VSCALE(n1, 1.0/(a2_len * a0_len));
   VSCALE(n2, 1.0/(a0_len * a1_len));

   THREE_VECTOR span;
   span.x = fabs(dot1(n0, box->bodyDiag));
   span.y = fabs(dot1(n1, box->bodyDiag));
   span.z = fabs(dot1(n2, box->bodyDiag));


   double minSpan = MIN(span.x, span.y);
   minSpan = MIN(minSpan, span.z);

   box->r2Inner = 0.25*minSpan*minSpan;

   // find needed shifts
   int intersection = 1;
   int shell = 0;
   box->lvSize = 0;
   while (box->itype != ORTHORHOMBIC && intersection == 1)
   {
      ++shell;
      int xShell = shell * (box->pbc & 1);
      int yShell = shell * (box->pbc & 2);
      int zShell = shell * (box->pbc & 4);
      intersection = 0;
      for (int ix=-xShell; ix<=xShell; ++ix)
         for (int iy=-yShell; iy<=yShell; ++iy)
            for (int iz=-zShell; iz<=zShell; ++iz)
            {
               if ( (abs(ix) != shell) && (abs(iy) != shell) && (abs(iz) != shell) )
                  continue;
               THREE_VECTOR ll = {ix, iy, iz};
               THREE_VECTOR rl = matrix_vector(box->h0, ll);
               double len = dot1(rl, rl);
               if (len <= 4*box->r2Outer) // could be more rigorous here
               {
                  intersection = 1;
                  if (box->lvCapacity < box->lvSize+1)
                  {
                     box->lvCapacity*=2;
                     box->lv = (THREE_INT*) ddcRealloc(box->lv, box->lvCapacity*sizeof(THREE_INT));
                  }
                  box->lv[box->lvSize].x = ix;
                  box->lv[box->lvSize].y = iy;
                  box->lv[box->lvSize].z = iz;
                  ++box->lvSize;
               }
            }
   }
}

/** Updates box parameters related to preduce and updates the preduce
 *  function pointers to select methods appropriate for the cell
 *  geometry.  This function should be called whenever the box geometry
 *  changes. */
void updatePreduce(BOX_STRUCT* box)
{
   findCellProperties(box);
   PsetMethod(box->pbc, &box->h0);
}

unsigned box_newAffineUpdateIndex(BOX_STRUCT* box)
{
   assert(box != NULL);

   unsigned index = box->nAffineUpdates;
   ++box->nAffineUpdates;

   box->affineUpdates = 
      ddcRealloc(box->affineUpdates, box->nAffineUpdates*sizeof(THREE_MATRIX));
   box->affineUpdates[index] = box->hinv;

   return index;
}

THREE_MATRIX box_getAffineUpdateHfac(BOX_STRUCT* box, unsigned index)
{
   if (box == NULL)
      box = _box_current;

   THREE_MATRIX hfac = box->h0;
   hfac = matrix_matrix(hfac, box->affineUpdates[index]);
   box->affineUpdates[index] = box->hinv;

   if (matrix_equal_tol(hfac, I_3x3, 1e-14))
      hfac=I_3x3;

   return hfac;
}



void misspelledKeywordError(void)
{
   if (getRank(0) != 0)
      return;

   printf("It appears that you are trying to use the bndcdn keyword in the\n"
         "BOX object, but you have misspelled it.  We have hit this problem\n"
         "so many times that we have changed the keyword name.  Please use\n"
         "pbc instead of bndcdn.\n");
   abortAll(-1);
}




/* Local Variables: */
/* tab-width: 3 */
/* End: */
