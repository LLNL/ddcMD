/* Specify a series of parallel planes; all normal to a box axis.
   These planes transect the simulation box.
   Axes must be Cartesian == box must be orthorhombic.
   N (>=2) planes delimit N disjoint regions.
   The plane positions are specified Before and After,
   and the delimited regions are transformed accordingly. */

#include "transectMorph.h"
#include <math.h>
#include <assert.h>
#include "simulate.h"
#include "transform.h"
#include "state.h"
#include "three_algebra.h" 
#include "object.h" 
#include "preduce.h"
#include "ddcMalloc.h"
#include "units.h"
#include "mpiUtils.h"

THREE_VECTOR transectNormal(TRANSFORM* transform, TRANSECTMORPH_PARMS * parms)
{
 // Return the 3-vector normal direction.
 // This is aligned with the index'th box axis.
 // DO NOT NORMALIZE the box axis length.
 SIMULATE* simulate = transform->parent;
 SYSTEM* sys = simulate->system;
 BOX_STRUCT *box = sys->box;

 THREE_VECTOR N;

 // which axis is normal aligned with?
 switch (parms->index)
        {
         case 0:
                 VSET(N,box->h0.xx,box->h0.yx,box->h0.zx);
                 break;
         case 1:
                 VSET(N,box->h0.xy,box->h0.yy,box->h0.zy);
                 break;
         case 2:
                 VSET(N,box->h0.xz,box->h0.yz,box->h0.zz);
                 break;
         default:
                 VSET(N,0.,0.,0.);
                 assert(parms->index>=0&&parms->index<=2);
        }
 
 return N;
}

////////////////////////////////////////////////////////////////////////
void transectMorph(TRANSFORM* transform)
{
 // Transform Before regions into After.
 // Just shift and stretch the space between adjacent planes.

 int taskNum = getRank(0);

   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;

   TRANSECTMORPH_PARMS* parms = transform->parms;

   THREE_VECTOR unitNormal = transectNormal(transform,parms);

   double halfL = 0.5*sqrt(VSQ(unitNormal));
   VSCALE(unitNormal,0.5/halfL);

 if (taskNum==0) printf(" %f %f %f box axis  ",unitNormal.x,unitNormal.y,unitNormal.z);
 if (taskNum==0) printf(" %f \n",halfL);
   
   unsigned nlocal = simulate->system->nlocal;


// which coordinate is getting morphed?
   double * ptr;
   switch (parms->index)
          {
           case 0: ptr = state->rx; break;
           case 1: ptr = state->ry; break;
           case 2: ptr = state->rz; break;
           default: assert(0==1);
          }
if (taskNum==0) printf(" here is halfL %f\n",halfL);
                  
   double tmpPos;
// double min=9e9, max=-9e9;
   int idid = 0;
   for (unsigned ii=0; ii<nlocal; ii++)
   {
      tmpPos = ptr[ii];
//    if (tmpPos<min) min=tmpPos;
//    if (tmpPos>max) max=tmpPos;

      if (ptr[ii]<parms->positionBefore[0])
         {
          double scaleB = (ptr[ii]-parms->positionBefore[0])
                         /(parms->positionBefore[0]+2.*halfL-parms->positionBefore[parms->nPlanes-1]);
          tmpPos = parms->positionAfter[0]
                   +scaleB
                        *(parms->positionAfter[0]+2.*halfL-parms->positionAfter[parms->nPlanes-1]);
          idid = 1;
         }

      for (int jj=1;jj<parms->nPlanes;jj++)
          {
//         printf(" %d %f %f\n",jj, ptr[ii],parms->positionBefore[jj]);
           if (ptr[ii]>=parms->positionBefore[jj-1] && ptr[ii]<parms->positionBefore[jj])
              {
               double scaleB = (ptr[ii]-parms->positionBefore[jj-1])
                              /(parms->positionBefore[jj]-parms->positionBefore[jj-1]);
               tmpPos = parms->positionAfter[jj-1]
                        +scaleB
                             *(parms->positionAfter[jj]-parms->positionAfter[jj-1]);
               idid = 2;
              }
           }
      if (ptr[ii]>parms->positionBefore[parms->nPlanes-1])
         {
          double scaleB = (ptr[ii]-parms->positionBefore[parms->nPlanes-1])
                         /(parms->positionBefore[0]+2.*halfL-parms->positionBefore[parms->nPlanes-1]);
          tmpPos = parms->positionAfter[parms->nPlanes-1]
                   +scaleB
                        *(parms->positionAfter[0]+2.*halfL-parms->positionAfter[parms->nPlanes-1]);
          idid = 3;
         }
      assert(idid>0);

      ptr[ii] = tmpPos;
   }

}

////////////////////////////////////////////////////////////////////////
void* transectMorph_parms(TRANSFORM* transform)
{
   TRANSECTMORPH_PARMS* parms = ddcMalloc(sizeof(TRANSECTMORPH_PARMS));
   int taskNum = getRank(0);

   OBJECT* obj = (OBJECT*) transform;

   transform->parms = parms;

   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system;
   BOX_STRUCT *box = sys->box;

   // Be sure box type is compatible with transectMorph transform
   assert (box->itype==ORTHORHOMBIC);

   // For now, require periodic boundary conditions in normal direction and all directions.
   assert (box->pbc==7); 

   // Specify index of normal: parallel to box axis # 0, 1, or 2
   object_get(obj, "index", &parms->index, INT, 1, "2");
   assert (parms->index>=0 && parms->index <=2);

   // Specify units of normal to be input: IN_UNITS(0) or IN_LATTICE(1)
   object_get(obj, "units", &parms->units, INT, 1, "0");
   assert (parms->units==0 || parms->units==1);

   // determine half the box length along normal
   THREE_VECTOR norm = transectNormal(transform,parms);
   double halfL = 0.5*sqrt(VSQ(norm));

   int ndumb;
   if (parms->units==IN_UNITS)
      { // here, positions are in length units
       parms->nPlanes = object_getv(obj, "positionBefore",  &parms->positionBefore,  DOUBLE, ABORT_IF_NOT_FOUND);
       ndumb          = object_getv(obj, "positionAfter",   &parms->positionAfter,   DOUBLE, ABORT_IF_NOT_FOUND);

      double length_convert = units_convert(1.0, NULL, "Angstrom");
      for (int i=0;i<parms->nPlanes;i++)
          {
           if (taskNum==0) printf("%d  %12.5e  %12.5e\n",i,parms->positionBefore[i],parms->positionAfter[i]);
           parms->positionBefore[i]/=length_convert;
           parms->positionAfter [i]/=length_convert;
          }
      }
   else if (parms->units==IN_LATTICE)
      { // here, positions are between -.5 and .5
       parms->nPlanes = object_getv(obj, "positionBefore",   parms->positionBefore,  DOUBLE, ABORT_IF_NOT_FOUND);
       ndumb          = object_getv(obj, "positionAfter",    parms->positionAfter,   DOUBLE, ABORT_IF_NOT_FOUND);
       for (int i=0;i<parms->nPlanes;i++)
           {
            parms->positionBefore[i] *= 2.*halfL;
            parms->positionAfter[i]  *= 2.*halfL;
           }
       assert(9<0); // disabled, not vetted
      }
   else
     { assert(3>3); }//nonexistent or erroneous option 

MPI_Barrier(COMM_LOCAL);

   // Verify same number of planes, Before and After.
   assert (ndumb==parms->nPlanes);

   // There must be at least two planes.
   assert (parms->nPlanes>=2);

   double origin;
   switch (parms->index)
          {
           case 0: origin = box->corner.x; break;
           case 1: origin = box->corner.y; break;
           case 2: origin = box->corner.z; break;
           default: assert(7==11);
          }

   // verify positional order of planes
   for (int i=1;i<parms->nPlanes;i++) assert(parms->positionBefore[i]>parms->positionBefore[i-1]);
   for (int i=1;i<parms->nPlanes;i++) assert(parms->positionAfter[i]>parms->positionAfter[i-1]);

   // verify that after planes don't 'pass' another on any wraparound
   if (parms->positionAfter[0]+2.*halfL<=parms->positionAfter[parms->nPlanes-1])
      if (getRank(0)==0) printf(" Transecting planes are crossing each other on wrapround?\n");
   if (!(parms->positionAfter[0]+2.*halfL > parms->positionAfter[parms->nPlanes-1]))
   {
      if (getRank(0) ==0) printf(" %e %e\n",parms->positionAfter[0]+2.*halfL , parms->positionAfter[parms->nPlanes-1]);
      assert(parms->positionAfter[0]+2.*halfL > parms->positionAfter[parms->nPlanes-1]);
   }

   // verify that Before planes lie inside unit cell.
   if (   parms->positionBefore[0]<origin
         || parms->positionBefore[parms->nPlanes-1]>=origin+2.*halfL)
   {
      printf(" transecting planes are not all inside the unit cell?\n");
      printf(" corner is at %f %f %f \n",box->corner.x,box->corner.y,box->corner.z);
   }
   assert (   parms->positionBefore[0]>=origin
         && parms->positionBefore[parms->nPlanes-1]<origin+2.*halfL);

   //cout << " *(#$(#$)(#$" << endl;

   return parms;

}
