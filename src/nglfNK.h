#ifndef NGLFNK_H
#define NGLFNK_H
#include "ddc.h"
#include "simulate.h"
#include "species.h"
typedef struct nglfMK_parms_st
{
   double P;  
   double T;
   double tau;
   THREE_VECTOR W;
} NGLFNK_PARMS;
NGLFNK_PARMS *nglfNK_parms(INTEGRATOR*integrator);
void nglfNK(DDC *ddc, SIMULATE *simulate, NGLFNK_PARMS *p);
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
