#include <stdio.h>
#include <stdlib.h>
#include "object.h"
#include "error.h"
#include "utilities.h"
#include "system.h"
#include "ddcMalloc.h"

typedef struct ZeroPotentialParms_st
{
   double rCut;
} ZERO_POTENTIAL_PARMS;



void *zeroPotential_parms(POTENTIAL *potential)
{
   ZERO_POTENTIAL_PARMS* parms = ddcMalloc(sizeof(ZERO_POTENTIAL_PARMS));
   OBJECT* obj = (OBJECT*) potential;

   object_get(obj, "rmax", &parms->rCut, WITH_UNITS, 1, "0.5", "Angstrom", NULL);

   return parms;
}
RCUT_TYPE *zeroPotentialCutoff(SYSTEM*sys, void* parmsV, int *n)
{
   ZERO_POTENTIAL_PARMS* parms = (ZERO_POTENTIAL_PARMS*) parmsV;
   int ncut = 1;
   static RCUT_TYPE rcut[2];
   rcut[0].value = parms->rCut;
   rcut[0].mode = RCUT_ALL;
   rcut[ncut] = rcut[0];
   *n = ncut;
   return rcut;
}
void zeroPotential(SYSTEM*sys, void *parms, ETYPE *e)
{
}
/* tab-width: 3 */
/* End: */
