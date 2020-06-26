#include "boxTransform.h"
#include <string.h>
#include "object.h"
#include "species.h"
#include "simulate.h"
#include "ddcMalloc.h"
#include "units.h"
#include "codata.h"
#include "box.h"

int getRank(int);


typedef struct boxTransform_parms_st
{
   THREE_MATRIX hNew;
} BOXTRANSFORM_PARMS;

void* boxTransform_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   BOXTRANSFORM_PARMS* parms = ddcMalloc(sizeof(BOXTRANSFORM_PARMS));

   object_get(obj, "hNew", &parms->hNew, WITH_UNITS, 9, "1 0 0 0 1 0 0 0 1", "l", NULL);

   return parms;
}


void boxTransform(TRANSFORM* transform)
{
   BOXTRANSFORM_PARMS* parms = transform->parms;
   SIMULATE* simulate = transform->parent;
   SYSTEM* sys = simulate->system;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;
   double* rx = state->rx;
   double* ry = state->ry;
   double* rz = state->rz;
   THREE_MATRIX hfac; 
   box_put(sys->box, HO, &parms->hNew);
   box_get(sys->box, HFAC, &hfac);
   if ( matrix_equal(hfac, I_3x3)) return ; 
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      THREE_VECTOR  rold={rx[ii],ry[ii],rz[ii]}; ;
      THREE_VECTOR r = matrix_vector(hfac,rold); 
      rx[ii]=r.x;
      ry[ii]=r.y;
      rz[ii]=r.z;
   }
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
