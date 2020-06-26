#ifndef ERFERFC_H
#define ERFERFC_H
#include <complex.h> 
typedef struct fdf_st { double f; double df;} FDF; 
typedef struct series_struct { double x0; double complex *c; int m;} SERIES_PARMS; 
double complex erfXOverXFast(double x, SERIES_PARMS *p);
double complex erfcFast(double x,SERIES_PARMS *p );
void erfcXOverX5(double *x2, double complex *fv,SERIES_PARMS *p);
SERIES_PARMS *erfcInit(int m,double x0);
SERIES_PARMS *erfInit(int m);
void freeSeriesParms(SERIES_PARMS *p);
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
