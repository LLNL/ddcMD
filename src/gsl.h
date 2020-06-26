#ifdef HAVE_GSL
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>
#else
#include <stddef.h>
#ifdef __cplusplus
extern "C" 
{
#endif
   typedef  double gsl_multifit_linear_workspace ; 
   typedef  struct { size_t size; double *v;}  gsl_vector ; 
   typedef  struct { size_t size_rows; size_t size_columns; double* m;}  gsl_matrix ; 
   double gsl_sf_fermi_dirac_mhalf(double x);
   double gsl_sf_fermi_dirac_half(double x);
   double gsl_sf_bessel_Kn(const int n, const double x);

   double gsl_vector_get (const gsl_vector * v, size_t i);
   gsl_matrix *gsl_matrix_alloc (size_t i, size_t j);
   gsl_vector *gsl_vector_alloc (size_t i);
   void gsl_matrix_set (gsl_matrix *m, size_t i, size_t j, double x);
   void gsl_vector_set (gsl_vector *v, size_t i, double x);
   void gsl_vector_free (gsl_vector * v);
   void gsl_matrix_free (gsl_matrix * v);
   gsl_multifit_linear_workspace *gsl_multifit_linear_alloc (size_t m,size_t n);
   void gsl_multifit_linear_free (gsl_multifit_linear_workspace *);
   int gsl_multifit_wlinear (const gsl_matrix * X, const gsl_vector * w, const gsl_vector * y, gsl_vector * c, gsl_matrix * cov, double * chisq, gsl_multifit_linear_workspace * work);
#ifdef __cplusplus
}
#endif
#endif


