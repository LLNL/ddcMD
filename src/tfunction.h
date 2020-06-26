#ifndef TFUNCTION_H
#define TFUNCTION_H

typedef struct tfunction_st
{
  int type;
  int rows;
  int cols;
  double x0;
  double dx;
  double *coef0;
  double *coef1;
  double *coef2;
} TFUNCTION;

TFUNCTION tfunc_init(const char* name);

void tfunc_f(TFUNCTION* this, double x, double* f1, double* f2);
void tfunc_df(TFUNCTION* this, double x, double* df1, double* df2);
void tfunc_fdf(TFUNCTION* this, double x, double* f, double* df);

double tfunc_xMax(TFUNCTION* this);
/** Returns the number of components in the value of the function.
 * Scalar functions have rank==1. */ 
unsigned tfunc_rank(TFUNCTION* this);


#endif // ifndef TFUNCTION_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
