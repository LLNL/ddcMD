#include <stdio.h>
#include <math.h>
typedef double (*PAIR_FUNCTION)(void* , void *parms, int, int, double,double *); 
static unsigned  npair_function = 0; 
static PAIR_FUNCTION pfunc_list[32] ; 
static void *pparms_list[32] ; 
void  pair_function_register( PAIR_FUNCTION pfunc, void *parms)
{	
	pfunc_list[npair_function] = pfunc; 
	pparms_list[npair_function] = parms; 
	npair_function++; 
}
double pair_function(void *sys, int i,int j,double r, double *dr)
{
	double p,dr_i; 
	p =*dr = 0.0; 
	for (unsigned k=0;k<npair_function;k++) 
	{
		p += pfunc_list[k](sys,pparms_list[k],i,j,r,&dr_i);
		*dr += dr_i;
	}
	return p; 
}
/*
{
	double v ; 
	v = sqrt(r); 
	*dr = 0.5/v; 
	return v; 
}
int main()
{
	double v,dr; 
	double parms = 2.0; 
	pair_function_register(myFunc,&parms);
	v= pair_function(0,0,4.0,&dr); 
	printf("%f %f\n",v,dr); 
}
*/


/* Local Variables: */
/* tab-width: 3 */
/* End: */
