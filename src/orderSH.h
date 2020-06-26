#ifndef ORDER_H
#define ORDER_H
#include "simulate.h"
#include "sph.h"
#include "functions.h"

enum ORDER_CLASS { SPHERICAL_HARMONIC, RHOK };

typedef struct orderSH_parms_st
{
   enum ORDER_CLASS itype;
   FUNCTION *function;
   int nL, L;
   SPH **sph;
   double Energy, dEnergydphi, phi;
   double Vo, Lo, r1o, r2o, r1, r2, lamda;
   int globalevalrate, localevalrate;
} ORDERSH_PARMS;

double*   orderGetQ(void);
int*      orderGetqN(void);
int*      orderGetC(void);
int*      orderGetLv(void);
int       orderGetnL(void);
double**  orderGetqnorm(void);
void writeqlocal(SIMULATE*simulate);


#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
