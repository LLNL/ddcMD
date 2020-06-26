#include "eam_fs.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "species.h"
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "ptiming.h"
#include "ddcMalloc.h"
#include "units.h"

typedef  struct {void *dummy;}  FS_EMBEDDING_PARMS; 

static double fs_embedding(FS_EMBEDDING_PARMS *parms, double rho, double *dv_dr);
static EP fs_pass1(FS_PASS_PARMS*parms, double r);
static EP fs_pass2(FS_PASS_PARMS*parms, double r);
static EP fs_series_pass1(FS_PASS_PARMS*parms, double r);
static EP fs_series_pass2(FS_PASS_PARMS*parms, double r);
void fs_pass0(FS_PASS_PARMS *parms, double r2, EP *ep,  EP *dp);


void eam_fs_parms(POTENTIAL *object, EAM_PARMS *parms)
{
	int i,j;
	FS_PASS_PARMS *pass_parms,**pp,*pi,*pj; 
	char **namelist,*mode; 

	int nspecies = parms->nspecies; 
	parms->pass_parms=ddcCalloc(nspecies*nspecies, sizeof(void*));
	parms->embedding_parms = ddcCalloc(nspecies, sizeof(void*));
	pp = (FS_PASS_PARMS **)parms->pass_parms; 
	for (i=0;i<nspecies;i++) 
	{
		parms->embedding_parms[i] = NULL;
		if (pp[i+i*nspecies]==NULL)
		   ddcMallocAligned((void*)&pp[i+i*nspecies],16,sizeof(FS_PASS_PARMS));
		for (j=i+1;j<nspecies;j++) 
		{
			if (pp[i+j*nspecies]==NULL)
			   ddcMallocAligned((void*)&pp[i+j*nspecies],16,sizeof(FS_PASS_PARMS));
			pp[j+i*nspecies ] = pp[i+j*nspecies];
		}
	}
	double length_convert = units_convert(1.0,"Angstrom",NULL); 
	double energy_convert = units_convert(1.0,"eV",NULL); 
	species_get(NULL,NAMELIST,(void *)&namelist);
	
	for (i=0;i<nspecies;i++) 
	{
		pass_parms = pp[i+i*nspecies];
		object_get((OBJECT *) object, namelist[i], &(pass_parms->a), DOUBLE, 6, "0.0"); 
		pass_parms->a *= energy_convert; 
		pass_parms->b *= energy_convert*energy_convert; 
		pass_parms->c *= length_convert; 
		pass_parms->l *= length_convert; 
		pass_parms->ro = 1.0*length_convert; 
		pass_parms->x = parms->rmax; 
	}
	for (i=0;i<nspecies;i++)
	for (j=i+1;j<nspecies;j++)
	{
		/* Combining rules */ 
		pi = pp[i+i*nspecies];
		pj = pp[j+j*nspecies];
		pass_parms = pp[i+j*nspecies];
		pass_parms->a = sqrt(pi->a*pj->a);
		pass_parms->b = sqrt(pi->b*pj->b);
		pass_parms->c = 0.25*(pi->c/pi->l + pj->c/pj->l)*(pi->l+pj->l);
		pass_parms->m = 0.5*(pi->m+pj->m);
		pass_parms->n = 0.5*(pi->n+pj->n);
		pass_parms->l = 0.5*(pi->l+pj->l);
		pass_parms->ro = 0.5*(pi->ro+pj->ro);
		pass_parms->x = 0.5*(pi->x+pj->x);
	}
	object_get((OBJECT *) object, "mode", &mode, STRING, 1, "ANALYTIC"); 
	if (strcmp(mode,"ANALYTIC") == 0 ) 
	{
		parms->embedding_function = (double (*)(void *, double , double *))fs_embedding; 
		parms->pass1_function = (EP (*)(void *,double))fs_pass1; 
		parms->pass2_function = (EP (*)(void *,double))fs_pass2; 
	}
	if (strcmp(mode,"SERIES") == 0  || strcmp(mode,"SERIES2") == 0 ) 
	{
		int k,l,cmax=0 ; 
		parms->embedding_function = (double (*)(void *, double, double *)) fs_embedding; 
		parms->pass1_function = (EP (*) (void *, double))fs_series_pass1; 
		parms->pass2_function = (EP (*) (void *, double))fs_series_pass2; 
        	double r_expansion; 
		object_get((OBJECT *) object, "r_expansion", &r_expansion, WITH_UNITS, 1, "0.0","Angstrom",NULL);
		for (i=0;i<nspecies;i++)
		{
			for (j=i;j<nspecies;j++)
			{
			   char keyword[1024]; 
			   double Cbody[64],Cdensity[64]; 
			   double DCbodyDr2[64],DCdensityDr2[64];
		
				for (l=0;l<64;l++) Cbody[l]=Cdensity[l]=DCbodyDr2[l]=DCdensityDr2[l]=0;
				sprintf(keyword,"%s-%s_2body",namelist[i],namelist[j]);
				if (object_testforkeyword((OBJECT*) object,keyword))
				{
					object_get((OBJECT *) object, keyword, Cbody, WITH_UNITS, 64, "0.0","eV",NULL); 
				}
				else
				{
					sprintf(keyword,"%s-%s_2body",namelist[j],namelist[i]);
					object_get((OBJECT *) object, keyword, Cbody, WITH_UNITS, 64, "0.0","eV",NULL); 
				}
				sprintf(keyword,"%s-%s_density",namelist[i],namelist[j]);
				if (object_testforkeyword((OBJECT*) object,keyword))
				{
					object_get((OBJECT *) object, keyword, Cdensity, WITH_UNITS,64,"0.0","eV^2",NULL);
				}
				else
				{
					sprintf(keyword,"%s-%s_density",namelist[j],namelist[i]);
					object_get((OBJECT *) object, keyword, Cdensity, WITH_UNITS,64,"0.0","eV^2",NULL);
				}
				if (strcmp(mode,"SERIES") == 0)
				{
					double alpha = units_convert(1.0,"Angstrom^-2",NULL); 
					double scale = 1.0; 
					cmax = 32; 
					for (l=0;l<cmax;l++) 
					{
						Cdensity[l] *= scale; 
						Cbody[l]    *= scale;  
						scale *= alpha;
					}
					for (l=0;l<cmax-1;l++) 
					{
						DCdensityDr2[l] = (l+1)*Cdensity[l+1]; 
						DCbodyDr2[l] =    (l+1)*Cbody[l+1] ;
					}
					DCdensityDr2[cmax-1]=0.0; 
					DCbodyDr2[cmax-1] =0.0  ;  
				}
				if (strcmp(mode,"SERIES2") == 0)
				{
					double scale_e,scale_d,alpha; 
					alpha = units_convert(1.0,"Angstrom^-2",NULL); 
					scale_e = units_convert(1.0,"Angstrom^8",NULL); 
					scale_d = units_convert(1.0,"Angstrom^4",NULL); 
					r_expansion = units_convert(3.1,"Angstrom",NULL); 
					double r2_expansion = r_expansion*r_expansion; 
					cmax = 35; 
					for (l=0;l<cmax;l++) 
					{
						
						Cdensity[l] *= scale_d; 
						Cbody[l]    *= scale_e;  
						scale_e *= alpha;
						scale_d *= alpha;
					}
					for (l=0;l<cmax-1;l++) 
					{
						DCdensityDr2[l] =  (l-2.0)*Cdensity[l] + r2_expansion*(l+1)*Cdensity[l+1]; 
						DCbodyDr2[l]    =  (l-4.0)*Cbody[l] + r2_expansion*(l+1)*Cbody[l+1]; 
					}
					DCdensityDr2[cmax-1]=DCbodyDr2[cmax-1] =0.0; 
				}
				pass_parms = pp[i+j*nspecies];
				pass_parms->r2_expansion = r_expansion*r_expansion; 
				if ( parms->SeriesOrder==PARTIAL_INTERLACE ) 
				{
					k = 2*cmax-1; 
    				for ( l=0;l<cmax;l++)
					{
						pass_parms->d[k-1] = -2.0*DCbodyDr2[l] ;    
						pass_parms->d[k  ] = -2.0*DCdensityDr2[l];  
						pass_parms->v[k-1] = Cbody[l];                          
						pass_parms->v[k  ] = Cdensity[l] ;                     
						k-=2; 
					}
				}
				if ( parms->SeriesOrder==FULL_INTERLACE ) 
				{
					k = 4*cmax-1; 
    				for ( l=0;l<cmax;l++)
					{
						pass_parms->v[k  ] = -2.0*DCdensityDr2[l];  
						pass_parms->v[k-1] = -2.0*DCbodyDr2[l] ;  
						pass_parms->v[k-2] = Cdensity[l] ;              
						pass_parms->v[k-3] = Cbody[l];                  
						k-=4; 
					}
				}
			}
		}
	}
}
double fs_embedding(FS_EMBEDDING_PARMS*parms, double rho, double *dv_dr)
{
	double v; 
	v = -sqrt(rho);
	*dv_dr=0.5/v;
	return v; 
}
void fs_pass0(FS_PASS_PARMS *parms, double r2, EP *ep,  EP *dp)
{
	double r,ri,dri,lr;
	ri = 1.0/sqrt(r2); 
	r = ri*r2; 
	dri =  1.0/(r-parms->x); 
	lr = log(r/parms->ro); 
	ep->e = parms->a*exp(parms->c*dri-parms->m*lr); 
	ep->p = parms->b*exp(parms->c*dri-parms->n*lr); 
	dp->e =-(parms->m*ri+parms->c*dri*dri)*ri*ep->e; 
	dp->p =-(parms->n*ri+parms->c*dri*dri)*ri*ep->p; 
}
EP fs_pass1(FS_PASS_PARMS *parms, double r2)
{
	EP ep; 
	double r,ri,dri,lr;
	ri = 1.0/sqrt(r2); 
	r = ri*r2; 
	dri =  1.0/(r-parms->x); 
	lr = log(r/parms->ro); 
	ep.e = parms->a*exp(parms->c*dri-parms->m*lr); 
	ep.p = parms->b*exp(parms->c*dri-parms->n*lr); 
	return ep; 
}
EP fs_pass2(FS_PASS_PARMS *parms, double r2)
{
	EP ep; 
	double r,lr,ri,dri;
	ri = 1.0/sqrt(r2); 
	r = ri*r2; 
	dri =  1.0/(r-parms->x); 
	lr = log(r/parms->ro); 
	ep.e = parms->a*exp(parms->c*dri-parms->m*lr); 
	ep.p = parms->b*exp(parms->c*dri-parms->n*lr); 
	ep.e *=-(parms->m/r+parms->c*dri*dri)*ri; 
	ep.p *=-(parms->n/r+parms->c*dri*dri)*ri; 
	return ep; 
}
EP fs_series_pass1(FS_PASS_PARMS *parms, double r2)
{
	EP ep; 
	double y,p0,t0,*d;
	int l; 
	d = parms->v; 
	y=parms->r2_expansion-r2; 
  	t0 = d[0];
  	p0 = d[1];
   	for (l=2;l<64;l+=2)
   	{
   		t0 = d[l]+y*t0;
   		p0 = d[l+1]  +y*p0;
   	}                                 
 	ep.p  = p0;
	ep.e  = t0;
	return ep; 
}
EP fs_series_pass2(FS_PASS_PARMS *parms, double r2)
{
	EP ep; 
	double y,p0,t0,*d;
	int l; 
	d = parms->d; 
	y=parms->r2_expansion-r2; 
  	t0 = d[0];
  	p0 = d[1];
   	for (l=2;l<64;l+=2)
   	{
   		t0 = d[l]+y*t0;
   		p0 = d[l+1]  +y*p0;
   	}                                 
 	ep.p  = p0;
	ep.e  = t0;
	return ep; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
