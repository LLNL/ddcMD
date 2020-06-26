#include <assert.h>
#include <math.h>
#define NMAX 32
static double _a[NMAX][NMAX];
double ErfcRootXOverRootX(double u,int n)
{
	double x = sqrt(u); 
	double uin,ui;
	ui = 1.0/u;
	uin = 1.0; 
	double f =0.0; 
	for (int k=1;k<=n;k++) 
	{
		uin *= ui; 
		f += _a[n][k]*uin;
	}
	f *= exp(-u)/sqrt(M_PI);
	f += _a[n][0]*erfc(x)*uin*sqrt(ui); 
	return f; 
}
void ErfcRootXOverRootXInit(int m)
{
	assert(m< NMAX);
	_a[0][0] =  1.0; 
	_a[0][1] = 0.0; 
	_a[1][0] = -0.5*_a[0][0];
	_a[1][1] = -1.0; 
	for (int n=1;n<m;n++)
	{
		_a[n+1][0] = -(2*n+1)/2.0*_a[n][0];
		_a[n+1][1] = -1.0*_a[n][1];
		for (int k=1;k<n;k++) _a[n+1][k+1]=(-k)*_a[n][k]-_a[n][k+1]; 
		_a[n+1][n+1] = (-n)*_a[n][n] - _a[n][0] ;
	}
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
