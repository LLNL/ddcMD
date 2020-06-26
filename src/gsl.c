#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "mpiUtils.h"
#ifndef HAVE_GSL
#include "gsl.h"
void  nogsl(void) 
{
   if (getRank(0)==0) printf("gsl library not available on  this platform"); fflush(stdout); 
   assert(1==0); 
}
double gsl_sf_fermi_dirac_mhalf(double x)
{
	nogsl(); 
	return 0.0; 
}
double gsl_sf_fermi_dirac_half(double x)
{
	nogsl(); 
	return 0.0; 
}
double gsl_sf_bessel_Kn(const int n, const double x)
{
        nogsl();
        return 0.0;
}
double gsl_vector_get (const gsl_vector * v, size_t i)
{
        nogsl();
        return 0.0;
}
gsl_matrix *gsl_matrix_alloc (size_t i, size_t j)
{
        nogsl();
        gsl_matrix *matrix=NULL; 
        return matrix;
}
gsl_vector *gsl_vector_alloc (size_t i)
{
        nogsl();
        gsl_vector *v=NULL;
        return v;
}
void gsl_matrix_set (gsl_matrix * m, size_t i, size_t j, double x)
{
        nogsl();
        return ;
}
void gsl_vector_set (gsl_vector * v, size_t i, double x)
{
        nogsl();
        return ;
}
void gsl_vector_free (gsl_vector * v)
{
        nogsl();
        return ;
}
void gsl_matrix_free (gsl_matrix * v)
{
        nogsl();
        return ;
}
gsl_multifit_linear_workspace *gsl_multifit_linear_alloc (size_t m,size_t n)
{
        nogsl();
	gsl_multifit_linear_workspace *workspace=NULL; 
        return  workspace; 
}
void gsl_multifit_linear_free (gsl_multifit_linear_workspace *work)
{
        nogsl();
        return; 
}
int gsl_multifit_wlinear (const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_vector * c, gsl_matrix * cov, double * chisq, gsl_multifit_linear_workspace * work)
{
        nogsl();
        return 1; 
}

#endif

int StopAppleComplainAboutNoSymbolsInThisFile=1;



