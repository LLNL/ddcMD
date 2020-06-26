#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#include "erfErfc.h"
#include "ddcMalloc.h"

#if defined(BGL) || defined(BGP) 
#define BG
#define ALIGN __align(16)
#include <mass.h> 
#define creal(x) __creal((x))
#define cimag(x) __cimag((x))
#else
#define ALIGN
#endif
#define NMAX  64
void freeSeriesParms(SERIES_PARMS *p)
{
	ddcFree(p->c); 
	ddcFree(p); 
}
double erfcSeries(int n,double *a, double x)
{
    double xn =1.0;
    double sum=0.0;
    for (int k=0;k<n;k++) { sum += a[k]*xn; xn *= x; }
    double fv = -2.0/sqrt(M_PI) * sum*exp(-x*x);
    return fv;
}
SERIES_PARMS *erfcInit(int m,double x0)   // Calculates the  coefficents for the m term talyor series approximations for erfc(x)  expanded about x0; 
{
    SERIES_PARMS *parms= ddcMalloc(sizeof(SERIES_PARMS)); 
    parms->c=NULL; 
    ddcMallocAligned((void*)&parms->c,16,m*sizeof(double complex));
    parms->m =m;  
    parms->x0 = x0; 
    double a[m][m];
    double b[m+1];
    b[0] = erfc(x0);
    a[1][0]=1.0; for (int k=1;k<m;k++) a[1][k]=0.0;
    b[1] = erfcSeries(1,a[1],x0);
    for (int n=2;n<m;n++)
    {
        for (int k=0;k<m;k++) a[n][k]=0.0;
        a[n][0] = a[n-1][1];
        for (int k=1;k<n;k++) a[n][k] = (k+1)*a[n-1][k+1]-2*a[n-1][k-1];
        b[n]=erfcSeries(n,a[n],x0);
    }
    b[m]=0.0;

    for (int n=0;n<m;n++)
    {
        double s = 1.0/tgamma(n+1.0);
        parms->c[m-n-1]=   s*(b[n]+I*b[n+1]);
    }
    return parms; 
}
double complex erfcFast(double x,SERIES_PARMS *p)
{
    double delta =  x-(p->x0);
    double complex   sum =  0.0;
    for (int k=0;k<p->m;k++)  sum = p->c[k]+delta*sum; //erfc(x)
    return sum;
}
void erfcXOverX4(double *x2,double complex *v,SERIES_PARMS *p)
{
	double complex f[4]; 
	double xinv[4],delta[4]; 
	double complex *c = p->c; 
	double m = p->m; 
	double x0 = p->x0; 
	
	for (int i=0;i<4;i++)
	{
		xinv[i] = 1.0/sqrt(x2[i]);
		delta[i] =  (x2[i]*xinv[i])-x0; 
	}
	f[0] = c[1]+delta[0]*c[0]; 
	f[1] = c[1]+delta[1]*c[0]; 
	f[2] = c[1]+delta[2]*c[0]; 
	f[3] = c[1]+delta[3]*c[0]; 
	for (int k=2;k<m;k++)  //erfc 
	{
	
		f[0] = c[k]+delta[0]*f[0]; 
		f[1] = c[k]+delta[1]*f[1]; 
		f[2] = c[k]+delta[2]*f[2]; 
		f[3] = c[k]+delta[3]*f[3]; 
	}
   	v[0] = (f[0]-I*creal(f[0])*xinv[0]*xinv[0])*xinv[0];
   	v[1] = (f[1]-I*creal(f[1])*xinv[1]*xinv[1])*xinv[1];
   	v[2] = (f[2]-I*creal(f[2])*xinv[2]*xinv[2])*xinv[2];
   	v[3] = (f[3]-I*creal(f[3])*xinv[3]*xinv[3])*xinv[3];
}
void erfcXOverX10(double *x2,double complex *v, SERIES_PARMS *p)
{
	double complex f[5]; 
	double xinv[10],delta[10]; 
	double complex *c = p->c; 
	double m = p->m; 
	double x0 = p->x0; 
	
	for (int i=0;i<5;i++)
	{
		xinv[i] = 1.0/sqrt(x2[i]);
		delta[i] =  (x2[i]*xinv[i])-x0; 
	}
	
	for (int i=0;i<2;i++)
	{
		f[0] = c[1]+delta[0]*c[0]; 
		f[1] = c[1]+delta[1]*c[0]; 
		f[2] = c[1]+delta[2]*c[0]; 
		f[3] = c[1]+delta[3]*c[0]; 
		f[4] = c[1]+delta[4]*c[0]; 
		for (int k=2;k<m;k++)  //erfc 
		{
			f[0] = c[k]+delta[0]*f[0]; 
			f[1] = c[k]+delta[1]*f[1]; 
			f[2] = c[k]+delta[2]*f[2]; 
			f[3] = c[k]+delta[3]*f[3]; 
			f[4] = c[k]+delta[4]*f[4]; 
		}
   	v[0] = (f[0]-I*creal(f[0])*xinv[0]*xinv[0])*xinv[0];
   	v[1] = (f[1]-I*creal(f[1])*xinv[1]*xinv[1])*xinv[1];
   	v[2] = (f[2]-I*creal(f[2])*xinv[2]*xinv[2])*xinv[2];
   	v[3] = (f[3]-I*creal(f[3])*xinv[3]*xinv[3])*xinv[3];
   	v[4] = (f[4]-I*creal(f[4])*xinv[4]*xinv[4])*xinv[4];
	}

}
void erfcXOverX5(double *x2,double complex *v, SERIES_PARMS *p)
{
	double complex f[5]; 
	double xinv[5],delta[5]; 
	double complex *c = p->c; 
	double m = p->m; 
	double x0 = p->x0; 
	
	for (int i=0;i<5;i++)
	{
		xinv[i] = 1.0/sqrt(x2[i]);
		delta[i] =  (x2[i]*xinv[i])-x0; 
	}
	f[0] = c[1]+delta[0]*c[0]; 
	f[1] = c[1]+delta[1]*c[0]; 
	f[2] = c[1]+delta[2]*c[0]; 
	f[3] = c[1]+delta[3]*c[0]; 
	f[4] = c[1]+delta[4]*c[0]; 
	for (int k=2;k<m;k++)  //erfc 
	{
	
		f[0] = c[k]+delta[0]*f[0]; 
		f[1] = c[k]+delta[1]*f[1]; 
		f[2] = c[k]+delta[2]*f[2]; 
		f[3] = c[k]+delta[3]*f[3]; 
		f[4] = c[k]+delta[4]*f[4]; 
	}
	((double *)v   ) [0]=creal(f[0])*xinv[0]; ((double *)v    )[1]=(creal(f[0])*xinv[0]-cimag(f[0]))*xinv[0]*xinv[0];
	((double *)(v+1))[0]=creal(f[1])*xinv[1]; ((double *)(v+1))[1]=(creal(f[1])*xinv[1]-cimag(f[1]))*xinv[1]*xinv[1];
	((double *)(v+2))[0]=creal(f[2])*xinv[2]; ((double *)(v+2))[1]=(creal(f[2])*xinv[2]-cimag(f[2]))*xinv[2]*xinv[2];
	((double *)(v+3))[0]=creal(f[3])*xinv[3]; ((double *)(v+3))[1]=(creal(f[3])*xinv[3]-cimag(f[3]))*xinv[3]*xinv[3];
	((double *)(v+4))[0]=creal(f[4])*xinv[4]; ((double *)(v+4))[1]=(creal(f[4])*xinv[4]-cimag(f[4]))*xinv[4]*xinv[4];
}
/*
 Function erfcXOverX returns a double complex value.  The real part is the function f(x2) = g(sqrt(x2))
  where g(x) = erfc(x)/x.   The imaginary piece is -g'(x)/x 
*/
double complex  erfcXOverX(double x2,SERIES_PARMS *p)
{
	double complex *c = p->c; 
	double m = p->m; 
	double x0 = p->x0; 
	double complex  v; 
	double xinv = 1.0/sqrt(x2);
	double delta =  x2*xinv-x0; 
	double complex   f =  c[0];  
	for (int k=1;k<m;k++)  //erfc 
	{
		f = c[k]+delta*f; 
	}
	((double *)&v)[0]=creal(f)*xinv; ((double *)&v)[1]=(creal(f)*xinv-cimag(f))*xinv*xinv;
//   	v = (f-I*creal(f)*xinv)*xinv;
	return v; 
}
SERIES_PARMS * erfInit(int m)   // Calculates the  coefficents for the m term talyor series approximations for erf(x^2)/x  expanded about 0; 
{
    SERIES_PARMS *parms= ddcMalloc(sizeof(SERIES_PARMS)); 
    parms->c=NULL; 
    ddcMallocAligned((void*)&parms->c,16,m*sizeof(double complex));
    parms->m = m;
    parms->x0 =0.0; 
    double b[m+1];
    double ai=M_2_SQRTPI; 
    b[0] = ai; 
    for (int i=1;i<m;i++)
    {
	ai = -ai*(2.0*i-1.0)/((2.0*i+1.0)); 
        b[i] = ai; 
    }
    b[m]=0.0;
    for (int n=0;n<m;n++)
    {
        double s = 1.0/tgamma(n+1.0);
        parms->c[m-n-1]=   s*(b[n]-2.0*I*b[n+1]);
    }
    return parms; 
}

double complex erfXOverXFast(double x,SERIES_PARMS *p)
{
    x = x*x; 
    double complex   sum =  0.0;
    for (int k=0;k<p->m;k++)  sum = p->c[k]+x*sum; //erf(x*x)/x
    return sum;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
