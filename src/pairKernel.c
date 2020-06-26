#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(BGL) || defined(BGP) 
#include <mass.h> 
#endif 
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "utilities.h"
#include "ptiming.h"
#include "opt.h"
extern int np,nspecies; 
extern double rcut,r2cut; 
extern double	sion_xx,sion_xy,sion_xz,sion_yy,sion_yz,sion_zz ;
extern FILE *ddcfile; 
int find_strips(OPT_BOX *box,OPT_STRIP *strip, int *nsi); 
void pairfcn(void *parms, int Tij, double r2, double *ep_ptr,double *dp_ptr)
{
	double ep,dp; 
	double y;
	y = r2cut-r2;
	ep = exp(-r2)*y; 
	dp = -2.0*ep - 2.0*exp(-r2); 
	*ep_ptr = ep; 
	*dp_ptr = dp; 
}

double block(void *parms , int Ti0, THREE_VECTOR ri, THREE_VECTOR *fi, double dFi, int first,int last,unsigned *type, THREE_VECTOR *rv, THREE_VECTOR *f,double *dF)
{
	double e2body=0.0;
	for (int j=first;j<=last;j++)
	{
		THREE_VECTOR r; 
		r.x = ri.x  - rv[j].x; 
		r.y = ri.y  - rv[j].y; 
		r.z = ri.z  - rv[j].z; 
		double r2  = r.x*r.x+r.y*r.y+r.z*r.z; 
		if (r2 < r2cut )  
		{
			unsigned Tij = Ti0 + type[j]; 
			double ep,dp; 
			pairfcn(parms,Tij,r2,&ep,&dp);
			double dFij = dFi+dF[j]; 
			e2body += dFij*ep;
			double fs =-(dFij*dp);
			double fxij=fs*r.x;
			double fyij=fs*r.y;
			double fzij=fs*r.z;
			fi->x +=  fxij;
			fi->y +=  fyij; 
			fi->z +=  fzij; 
			f[j].x -=  fxij; 
			f[j].y -=  fyij; 
			f[j].z -=  fzij;
/*
			sion_xx -= fxij*r.x; 
			sion_xy -= fxij*r.y; 
			sion_xz -= fxij*r.z; 
			sion_yy -= fyij*r.y; 
			sion_yz -= fyij*r.z; 
			sion_zz -= fzij*r.z; 
*/
		}
	}
	return e2body; 
}
double  pairKernel(void *parms, OPT_BOX *box, OPT_PARTICLE *rv, THREE_VECTOR *simage, unsigned *type, double *dF, THREE_VECTOR *f)
{
	THREE_VECTOR ri,s; 
	double e2body;
	int i,k;
	int image_start,ns,first_i,last_i,first_j,last_j; 
	OPT_STRIP strip[27]; 
	ns= find_strips(box,strip,&image_start);
	e2body=0.0; 
	s = *simage;
	first_i = box->first; 
	last_i = box->last; 
	for (i=first_i;i<=last_i;i++)
	{
		double dFi = dF[i]; 
		THREE_VECTOR fi = {0.0,0.0,0.0};
		ri.x = rv[i].x; ri.y = rv[i].y; ri.z = rv[i].z;
		int Ti0 = type[i]*nspecies; 
		e2body += block(parms, Ti0,ri, &fi, dFi, i+1,last_i,type, (THREE_VECTOR *)rv, f,dF);
		for (k =0;k<image_start;k++) 
		{
			first_j = strip[k].first; 
			last_j = strip[k].last; 
			e2body += block(parms, Ti0, ri, &fi, dFi, first_j,last_j,type, (THREE_VECTOR *)rv, f,dF);
		}
		ri.x += s.x ; ri.y += s.y ; ri.z += s.z ; 
		for (k=image_start;k<ns;k++) 
		{
			first_j = strip[k].first; 
			last_j = strip[k].last; 
			e2body += block(parms, Ti0, ri, &fi, dFi, first_j,last_j,type, (THREE_VECTOR *)rv, f,dF);
		}
		f[i].x +=  fi.x; 
		f[i].y +=  fi.y; 
		f[i].z +=  fi.z;
	}
	return e2body; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
