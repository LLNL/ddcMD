#include "eam.h" 
#include "eam_fs.h" 
#include "eam_exp.h"
#include "eam_sc.h"
#include "eam_at.h"
#include "eam_tabular.h"
#include "eam_rational.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "ptiming.h"
#include "species.h"
#include "system.h"
#include "ddcMalloc.h"
#include "units.h"
#include "expandbuffer.h"
#include "hpmWrapper.h"

/* EAM potential routines */
EAM_PARMS *eam_parms(POTENTIAL *object)
{
	char **namelist; 
	int nspecies=system_getNspecies(NULL);
	species_get(NULL,NAMELIST,(void *)&namelist);

	EAM_PARMS* parms = ddcMalloc(sizeof(EAM_PARMS));
	parms->nspecies=nspecies; 
	char *form; 
	object_get((OBJECT *) object, "form", &form, STRING, 1, "exp");
	object_get((OBJECT *) object, "rmax", &parms->rmax, WITH_UNITS, 1, "0.0","Angstrom",NULL);
	parms->pass_parms=NULL;
	parms->embedding_parms=NULL;
	parms->SeriesOrder = PARTIAL_INTERLACE; 
	if (strcmp(form,"EXP") == 0)      eam_exp_parms(object,parms); 
	if (strcmp(form,"FS") == 0)       eam_fs_parms(object,parms);
	if (strcmp(form, "SC") == 0)      eam_sc_parms(object, parms);
	if (strcmp(form, "AT") == 0)      eam_at_parms(object, parms);
	if (strcmp(form, "TABULAR") == 0) eam_tabular_parms(object, parms);
	if (strcmp(form, "RATIONAL") == 0) eam_rational_parms(object, parms);
	
	return parms;
}
RCUT_TYPE *eamCutoff(SYSTEM*sys, EAM_PARMS*parms, int *n)
{
	static RCUT_TYPE rcut[2];
	static int ncut = 1;
	rcut[0].value = parms->rmax;
	rcut[0].mode = RCUT_ALL;
	rcut[ncut] = rcut[0];
	*n = ncut;
	return rcut;
}
void eam(SYSTEM*sys, EAM_PARMS*parms, ETYPE* eStuff)
{
	/* EAM vars */
 	static double *rho = NULL, *Fi = NULL, *dFi_drho = NULL; 
	HPM_Start("eam");
	profile(P_KER, START);

	double* fx = sys->collection->state->fx;
	double* fy = sys->collection->state->fy;
	double* fz = sys->collection->state->fz;
	int*  type = sys->collection->state->atomtype;
	double* energy = sys->collection->state->potentialEnergy;
	double* virial = sys->collection->state->virial;
	THREE_SMATRIX *sion = sys->collection->state->sion;

	NBR* nbr = sys->neighbor;
	int nlocal = sys->nlocal;
	int nion = sys->nion; 
	NPARTICLE* particles = sys->neighbor->particles;
	rho = (double *)ExpandBuffers((void *)rho, sizeof(double), nion, 512, LOCATION("eam"), "rho");
	Fi = (double *)ExpandBuffers((void *)Fi, sizeof(double), nlocal, 512, LOCATION("eam"), "Fi");
	dFi_drho = (double *)ExpandBuffers((void *)dFi_drho, sizeof(double), nion, 512, LOCATION("eam"), "dFi_drho");
	int IndexRmax = -1;
	for (int i = 0; i < nbr->nc; i++) if (fabs(nbr->rcut[i].value - parms->rmax) < 1e-8) IndexRmax = i;
	if (IndexRmax == -1 )
	{
		printf("Error in eam. Index rcut\n");
		exit(1);
	}

	for (int i = 0; i < nlocal; i++) rho[i] = Fi[i] = dFi_drho[i] = 0.0;

	/* compute densities for all local atoms */
	// unused double rcut = parms->rmax;  
	// unused	int nc = 4*parms->cmax; 
	int nspecies = parms->nspecies; 
	double engy2 =0.0;
	double engyFi = 0.0;
	profile(EAM_PASS1, START);
	for (int i = 0; i < nlocal; i++) 
	{
	   double energy_i = 0;
	   double rho_i=0.0; 
	   PAIRS* pij = particles[i].ifirst[IndexRmax];
	   int Ti = type[i]&0xffff; 
	   while (pij != NULL )
	   {
	      EP ep; 
	      double r2; 
	      int j = pij->j; 
	      RPAIR* p = pij->p;
	      int Tj = type[j]&0xffff; 
	      r2 = p->r*p->r; 
	      ep = parms->pass1_function(parms->pass_parms[Ti+Tj*nspecies],r2); 
	      if (j<nlocal)
	      {
				engy2 += ep.e; 
				energy[j] += 0.5*ep.e;
				energy_i += 0.5*ep.e;
	      }
	      rho_i += ep.p;
	      rho[j] += ep.p;
	      pij = pij->ilink;
	   }	
	   rho[i] += rho_i; 
	   energy[i] += energy_i;
	}
	for (int i = nlocal; i < nion; i++) 
	{
	   double rho_i=0.0;
	   double energy_i=0;
	   PAIRS* pij = particles[i].ifirst[IndexRmax];
	   int Ti = type[i]&0xffff; 
	   while (pij != NULL )
	   {
	      int j = pij->j; 
	      if (j < nlocal)   
	      {
				RPAIR* p = pij->p;
				int Tj = type[j]&0xffff; 
				double r2 = p->r*p->r; 
				EP ep = parms->pass1_function(parms->pass_parms[Ti+Tj*nspecies],r2); 
				engy2 += ep.e ; 
				energy[j] += 0.5*ep.e;
				energy_i += 0.5*ep.e;
				rho_i += ep.p;
				rho[j] += ep.p;
	      }
	      pij = pij->ilink;
	   }	
	   rho[i] += rho_i; 
	   energy[i] += energy_i;
	}
	profile(EAM_PASS1, END);
	
	/* compute  F, dF/drho for each local and remote atom */
	/* compute F(rho) and dF/drho  arrays   */
	profile(EAM_EMBED, START);
	for (int i = 0; i < nlocal; i++)  
	{
		int Ti = type[i]&0xffff; 
		double Fi = parms->embedding_function(parms->embedding_parms[Ti], rho[i], dFi_drho+i);
		energy[i] += Fi; 
		engyFi += Fi; 
	}
	profile(EAM_EMBED, END);
	for (int i = nlocal; i < nion; i++) dFi_drho[i]=0.0; 
	/* compute pair potential */
	profile(EAM_PASS2, START);
	for (int i = 0; i < nlocal; i++)
	{
		EP ep; 
		double p2;
		double virial_i,fx_i,fy_i,fz_i;
		fx_i=fy_i=fz_i=0.0; 
		virial_i=0.0; 
		THREE_SMATRIX sion_i = {0.0,0.0,0.0,0.0,0.0,0.0}; 
		PAIRS* pij = particles[i].ifirst[IndexRmax];
		int Ti = type[i]&0xffff; 
		while (pij  != NULL)
		{
			int j = pij->j;
			RPAIR* p = pij->p;
			int Tj = type[j]&0xffff; 
			ep = parms->pass2_function(parms->pass_parms[Ti+Tj*nspecies],p->r*p->r); 
			p2 =-(ep.e + ep.p*(dFi_drho[i]+dFi_drho[j]));
			{
				double x,y,z,fpx, fpy, fpz;
				x = p->x; y=p->y; z=p->z; 
				fpx = p2*x;
				fpy = p2*y;
				fpz = p2*z;
			
				fx_i += fpx;
				fy_i += fpy;
				fz_i += fpz;
				fx[j] -= fpx;
				fy[j] -= fpy;
				fz[j] -= fpz;
				virial_i += 0.1666666666666667*(fpx*x+fpy*y+fpz*z);
				virial[j] += 0.1666666666666667*(fpx*x+fpy*y+fpz*z);
				double  fpxx = fpx*x;
				double  fpxy = fpx*y;
				double  fpxz = fpx*z;
				double  fpyy = fpy*y;
				double  fpyz = fpy*z;
				double  fpzz = fpz*z;

				sion_i.xx -= 0.5*fpxx;
				sion_i.xy -= 0.5*fpxy;
				sion_i.xz -= 0.5*fpxz;
				sion_i.yy -= 0.5*fpyy;
				sion_i.yz -= 0.5*fpyz;
				sion_i.zz -= 0.5*fpzz;
			
				sion[j].xx -= 0.5*fpxx;
				sion[j].xy -= 0.5*fpxy;
				sion[j].xz -= 0.5*fpxz;
				sion[j].yy -= 0.5*fpyy;
				sion[j].yz -= 0.5*fpyz;
				sion[j].zz -= 0.5*fpzz;

				eStuff->virial.xx += fpxx;
				eStuff->virial.xy += fpxy;
				eStuff->virial.xz += fpxz;
				eStuff->virial.yy += fpyy;
				eStuff->virial.yz += fpyz;
				eStuff->virial.zz += fpzz;
			}
			pij = pij->ilink;
		}
		fx[i] += fx_i; 
		fy[i] += fy_i; 
		fz[i] += fz_i; 

		virial[i] += virial_i; 

		sion[i].xx += sion_i.xx;
		sion[i].xy += sion_i.xy;
		sion[i].xz += sion_i.xz;
		sion[i].yy += sion_i.yy;
		sion[i].yz += sion_i.yz;
		sion[i].zz += sion_i.zz;
	}
	for (int i = nlocal; i < nion; i++)
	{	
		EP ep; 
		double p2;
		double virial_i,fx_i,fy_i,fz_i;
		fx_i=fy_i=fz_i=0.0;
		virial_i=0.0; 
		THREE_SMATRIX sion_i = {0.0,0.0,0.0,0.0,0.0,0.0}; 
		int Ti = type[i]&0xffff; 
		PAIRS* pij = particles[i].ifirst[IndexRmax];
		while (pij  != NULL)
		{
			int j = pij->j;
			if (j < nlocal)
			{
				RPAIR* p = pij->p;
				int Tj = type[j]&0xffff; 
				ep = parms->pass2_function(parms->pass_parms[Ti+Tj*nspecies],p->r*p->r); 
				p2 =-(ep.p*(dFi_drho[j]));
				//double r2=p->r*p->r; 
				{
				double x,y,z,fpx, fpy, fpz;
				x = p->x; y=p->y; z=p->z; 
				fpx = p2*x;
				fpy = p2*y;
				fpz = p2*z;
			
				fx_i += fpx;
				fy_i += fpy;
				fz_i += fpz;
				fx[j] -= fpx;
				fy[j] -= fpy;
				fz[j] -= fpz;
				virial_i += 0.1666666666666667*(fpx*x+fpy*y+fpz*z);
				virial[j] += 0.1666666666666667*(fpx*x+fpy*y+fpz*z);
				double  fpxx = fpx*x;
				double  fpxy = fpx*y;
				double  fpxz = fpx*z;
				double  fpyy = fpy*y;
				double  fpyz = fpy*z;
				double  fpzz = fpz*z;
			
				sion_i.xx -= 0.5*fpxx;
				sion_i.xy -= 0.5*fpxy;
				sion_i.xz -= 0.5*fpxz;
				sion_i.yy -= 0.5*fpyy;
				sion_i.yz -= 0.5*fpyz;
				sion_i.zz -= 0.5*fpzz;
			
				sion[j].xx -= 0.5*fpxx;
				sion[j].xy -= 0.5*fpxy;
				sion[j].xz -= 0.5*fpxz;
				sion[j].yy -= 0.5*fpyy;
				sion[j].yz -= 0.5*fpyz;
				sion[j].zz -= 0.5*fpzz;

				eStuff->virial.xx += fpxx;
				eStuff->virial.xy += fpxy;
				eStuff->virial.xz += fpxz;
				eStuff->virial.yy += fpyy;
				eStuff->virial.yz += fpyz;
				eStuff->virial.zz += fpzz;
				}
			}
			pij = pij->ilink;
		}
		fx[i] += fx_i; 
		fy[i] += fy_i; 
		fz[i] += fz_i; 
		virial[i] += virial_i; 
		sion[i].xx += sion_i.xx;
		sion[i].xy += sion_i.xy;
		sion[i].xz += sion_i.xz;
		sion[i].yy += sion_i.yy;
		sion[i].yz += sion_i.yz;
		sion[i].zz += sion_i.zz;
	}
	profile(EAM_PASS2, END);
//	printf("%d: e=%f %f\n",getRank(0),engy2,engyFi);
	eStuff->eion += engyFi + engy2;
	profile(P_KER, END);
	HPM_Stop("eam");
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
