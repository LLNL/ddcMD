#ifndef EAM_RATIONAL__
#define EAM_RATIONAL__

#include "potential.h"
#include "eam.h"

typedef struct {
  int degree;
  double *coeff,(*val_and_deriv_coeff)[2];
} polynomial;
typedef struct {
  double cutoff;
  polynomial numerator,denominator;
} ratexpr;
typedef struct {
  ratexpr rho_fun,phi_fun;
} RATIONAL_PASS_PARMS;
typedef struct {
  ratexpr F_fun;
} RATIONAL_EMBEDDING_PARMS;


void eam_rational_parms(POTENTIAL *object, EAM_PARMS *parms);

double eval_rational(int num_deg,const double num_coeff[],
		     int den_deg,const double den_coeff[],
		     double x,double deriv[1]);


#endif
