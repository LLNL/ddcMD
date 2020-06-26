#ifndef PREDUCE_H
#define PREDUCE_H

#include "three_algebra.h"

#ifdef __cplusplus
extern "C" 
{
#endif

extern void (*nearestImage)(double*, double*, double*);
extern void (*backInBox)(double*, double*, double*);
extern void (*nearestImage_fast)(double*, double*, double*);
extern void (*backInBox_fast)(double*, double*, double*);

void Pset(double *hptr, double *hiptr, int *bndptr, double* r2Safe, int* nShifts, THREE_INT** lShifts);
void PsetMethod(unsigned pbc, THREE_MATRIX* h);

#ifdef __cplusplus
}
#endif

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
