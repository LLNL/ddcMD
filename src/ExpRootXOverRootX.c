#include <math.h>
#include <assert.h>
#include <stdio.h>
#define NMAX 32
static double _a[NMAX][NMAX];
double ExpRootXOverRootX(double u,int n)
{
	double x = sqrt(u); 
	double xn =1.0; 
	double s=-2.0; 
	double f =0.0; 
	for (int k=0;k<=n;k++) 
	{
		f += _a[n][k]*xn;
		s *= -0.5; 
		xn *= x; 
	}
	f *= s*exp(-x)/pow(x,2*n+1);
	return f; 
}
void ExpRootXOverRootXInit(int m)
{
	assert(m< NMAX);
	for (int n=0;n<m;n++) for (int k=0;k<m;k++) _a[n][k] = 0.0; 
	_a[0][0] =  1.0; 
	for (int n=0;n<m;n++)
	{
		_a[n+1][0] = (2*n+1)*_a[n][0];
		for (int k=1;k<=n;k++) _a[n+1][k]=(2*n+1-k)*_a[n][k]+_a[n][k-1]; 
		_a[n+1][n+1]=_a[n][n]; 
	}
//	for (int n=0;n<m;n++) {for (int k=0;k<=n;k++) printf("%f ",_a[n][k]); printf("\n"); }
	
}
/*
int main()
{
	ExpRootXOverRootXInit(6);
	for (double x=0.001;x<2.0;x+=0.001)
	{
		double rx = sqrt(x); 
		printf("%f %e %e %e\n",x,exp(-rx)/rx,ExpRootXOverRootX(x,1),ExpRootXOverRootX(x,2)); 
	}
}*/


/* Local Variables: */
/* tab-width: 3 */
/* End: */
