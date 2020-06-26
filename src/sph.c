#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "three_algebra.h"
//#include "error.h"
#include  "sph.h"
#include "ddcMalloc.h"
typedef struct complex_st
{
	double real, imag;
} COMPLEX;

typedef struct ylm_st
{
	COMPLEX value, x, y, z;
} YLM;
double cut=1e-8; 

#define LMAX  16
void sph(SPH*parms, double x, double y, double z, YLM*ylm);
SPH *sphinit(int L)
{
	int m;
	SPH *parms;
	parms = (SPH *) ddcMalloc(sizeof(SPH));
	parms->L = L;
	//double A = parms->A = sqrt((2.0*L + 1.0)/(4.0*M_PI));
	double *B = parms->B = (double *)ddcMalloc((L + 1)*sizeof(double));
	double *C = parms->C = (double *)ddcMalloc((L + 1)*sizeof(double));
	for (m = 0; m <= L; m++)
	{
		B[m] = 1.0/sqrt((L + 1.0 - m)*(L + m));
		C[m] = -1.0/((L - m)*(L + m + 1.0));
	}
	return parms;
}

void sph(SPH*parms, double x, double y, double z, YLM*ylm)
{
	double A,*B, *C;
	double u, v, rxy, rxyi, ri, z1, xz, yz, zz, p[LMAX + 1];
	double a, c, dpdcos;
	int l, m, L;
	L = parms->L;
	A = parms->A;
	B = parms->B;
	C = parms->C;
	rxy = sqrt(x*x + y*y);
	ri = 1.0/sqrt(x*x + y*y + z*z);
	if (rxy*ri < cut)
	{
		x *= ri;
		y *= ri;
		z *= ri;
		c =-0.5*A*L*(L+1.0)*ri;
		ylm[0].value.real = A; 
		ylm[0].value.imag = 0.0;  
		ylm[0].x.real= c*x;
		ylm[0].x.imag = 0.0;  
		ylm[0].y.real= c*y;
		ylm[0].y.imag = 0.0;  
		ylm[0].z.real= 0.0;
		ylm[0].z.imag = 0.0;  
		if (L==0) return; 
		c =0.5*A*sqrt(L*(L+1.0));
		ylm[1].value.real =-c*x; 
		ylm[1].value.imag =-c*y;  
		ylm[1].x.real= -c*ri;
		ylm[1].x.imag = 0.0;  
		ylm[1].y.real = 0.0;  
		ylm[1].y.imag= -c*ri;
		ylm[1].z.real= c*ri*x;
		ylm[1].z.imag =c*ri*y;
		if (L==1) return; 
		c =0.25*A*sqrt((L-1.0)*L*(L+1.0)*(L+2.0))*ri;
		ylm[2].value.real = 0.0; 
		ylm[2].value.imag = 0.0;  
		ylm[2].x.real= c*x;
		ylm[2].x.imag = c*y;  
		ylm[2].y.real =-c*y;  
		ylm[2].y.imag=  c*x;
		ylm[2].z.real= 0.0;
		ylm[2].z.imag =0.0; 
		if (L==2) return; 
		for (l=2;l<=L;l++) 
		{
			ylm[l].value.real = 0.0; 
			ylm[l].value.imag = 0.0;  
			ylm[l].x.real= 0.0;
			ylm[l].x.imag =0.0; 
			ylm[l].y.real= 0.0;
			ylm[l].y.imag =0.0; 
			ylm[l].z.real= 0.0;
			ylm[l].z.imag =0.0; 
		}
		return; 
	}
	else
	{
		rxyi = 1.0/rxy;
	}
	u = x*rxyi;
	v = y*rxyi;
	x *= ri;
	y *= ri;
	z *= ri;
	z1 = rxy*ri;
	a = z/z1;
	ylm[0].value.real = A;
	ylm[0].value.imag = 0.0;
	p[L + 1] = 0.0;
	p[L] = 1.0;
	for (m = 1; m <= L; m++)
	{
		ylm[m].value.real = (u*ylm[m - 1].value.real - v*ylm[m - 1].value.imag)*B[m];
		ylm[m].value.imag = (v*ylm[m - 1].value.real + u*ylm[m - 1].value.imag)*B[m];
		p[L] *= -z1*(2.0*m - 1.0);
	}
	p[L - 1] = -a*p[L];
	for (m = L - 1; m > 0; m--)
	{
		p[m - 1] = (m*(2.0*a)*p[m] + p[m + 1])*C[m - 1];
	}
	xz = -x*a*ri;
	yz = -y*a*ri;
	zz = z1*ri;
	for (m = 0; m <= L; m++)
	{
		dpdcos = -m*a*p[m] - p[m + 1];
		c = dpdcos*xz;
		ylm[m].x.real = c*ylm[m].value.real;
		ylm[m].x.imag = c*ylm[m].value.imag;
		c = dpdcos*yz;
		ylm[m].y.real = c*ylm[m].value.real;
		ylm[m].y.imag = c*ylm[m].value.imag;
		c = dpdcos*zz;
		ylm[m].z.real = c*ylm[m].value.real;
		ylm[m].z.imag = c*ylm[m].value.imag;

		ylm[m].value.real *= p[m];
		ylm[m].value.imag *= p[m];

		ylm[m].x.real += m*(v*rxyi)*(ylm[m].value.imag);
		ylm[m].x.imag -= m*(v*rxyi)*(ylm[m].value.real);
		ylm[m].y.real -= m*(u*rxyi)*(ylm[m].value.imag);
		ylm[m].y.imag += m*(u*rxyi)*(ylm[m].value.real);
	}
}
/*
int main()
{
	int i,l; 
	double r,x,y,z,theta,phi;
	SPH *parm; 
	YLM s[7]; 
	scanf("%d %lf %lf %lf %le",&l,&r,&theta,&phi,&cut);
	for (i=0;i<15;i++)
	{
	x = r*sin(theta)*cos(phi);
	y = r*sin(theta)*sin(phi);
	z = r*cos(theta);
	parm = sphinit(l);
	sph(parm,x,y,z,s);
	printf("%14.12f %14.12f %14.12f %14.12f %14.12f %14.12f %14.12f %14.12f %14.12f\n",theta,s[2].value.real,s[2].value.imag, s[2].x.real,s[2].x.imag,s[2].y.real,s[2].y.imag,s[2].z.real,s[2].z.imag);
	theta *= 0.5; 
	}
	
}
*/


/* Local Variables: */
/* tab-width: 3 */
/* End: */
