//To do:
//
// 1.  Calculation of embedding parameters isn't species dependent (this
// may be ok since this may be an alternate implementation of
// Finnis-Sinclair form).

#include "eam_exp.h" 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "utilities.h"
#include "ptiming.h"
#include "eam.h"
#include "potential.h"
#include "ddcMalloc.h"

typedef  struct exp_embedding_parms { double rho_e, E_c, alpha, beta, gamma, rho_e_inv, ab, gb;} EXP_EMBEDDING_PARMS ;
typedef  struct exp_pass_parms { double v[64],d[64],r_expansion,r_e,beta,gamma,f_e,phi_e,r_e_inv;} EXP_PASS_PARMS ;

static double exp_embedding(EXP_EMBEDDING_PARMS*parms, double rho, double *dv_dr);
static EP exp_pass1(EXP_PASS_PARMS*parms, double r2);
static EP exp_pass2(EXP_PASS_PARMS*parms, double r2);

void eam_exp_parms(POTENTIAL *object, EAM_PARMS *parms)
{
	int i,j,nspecies;
	double atomvolume,r_e,rho_e,phi_e,alpha,beta,gamma,E_c; 
	EXP_PASS_PARMS *pass_parms,**pp; 

	nspecies = parms->nspecies; 
	parms->embedding_parms = ddcCalloc(nspecies, sizeof(EXP_EMBEDDING_PARMS*));
	parms->pass_parms = ddcCalloc(nspecies*nspecies,sizeof(EXP_PASS_PARMS*));
	object_get((OBJECT *) object, "atomvolume", &atomvolume, WITH_UNITS, 1, "1.0","Angstrom^3",NULL);
	object_get((OBJECT *) object, "phi_e", &phi_e, WITH_UNITS, 1, "0.0","eV",NULL);
	object_get((OBJECT *) object, "r_e", &r_e, WITH_UNITS, 1, "0.0","Angstrom",NULL);
	object_get((OBJECT *) object, "rho_e", &rho_e, DOUBLE, 1, "0.0");
	object_get((OBJECT *) object, "alpha", &alpha, DOUBLE, 1, "0.0");
	object_get((OBJECT *) object, "beta", &beta, DOUBLE, 1, "0.0");
	object_get((OBJECT *) object, "gamma", &gamma, DOUBLE, 1, "0.0");
	object_get((OBJECT *) object, "E_c", &E_c, WITH_UNITS, 1, "0.0","eV",NULL);
	EXP_EMBEDDING_PARMS embedding_parms; 
	embedding_parms.rho_e = rho_e=E_c/atomvolume;
	embedding_parms.rho_e_inv = 1.0/embedding_parms.rho_e;  
	embedding_parms.E_c = E_c;
	embedding_parms.alpha = alpha;
	embedding_parms.beta = beta;
	embedding_parms.gamma = gamma;
	embedding_parms.ab = alpha/beta;
	embedding_parms.gb = gamma/beta;
	for (i=0;i<nspecies;i++)
	{
		parms->embedding_parms[i] = ddcCalloc(1, sizeof(EXP_EMBEDDING_PARMS));
		*((EXP_EMBEDDING_PARMS *)(parms->embedding_parms[i]))=embedding_parms; 
	}
	pp = (EXP_PASS_PARMS **)parms->pass_parms; 
	for (i=0;i<nspecies;i++) 
	{
		if (pp[i+i*nspecies]==NULL) pp[i+i*nspecies] = ddcMalloc(sizeof(EXP_PASS_PARMS));
		for (j=i+1;j<nspecies;j++) 
		{
			if (pp[i+j*nspecies] == NULL ) pp[j+nspecies*i]= pp[i+j*nspecies]=ddcMalloc(sizeof(EXP_PASS_PARMS));
		}
	}
	for (i=0;i<nspecies;i++)
	{
		pass_parms=pp[i+i*nspecies]; 
		pass_parms->r_e = r_e; 
		pass_parms->r_e_inv = 1.0/r_e; 
		pass_parms->beta = beta; 
		pass_parms->gamma = gamma; 
		pass_parms->f_e = rho_e/12.0;
		pass_parms->phi_e = E_c/6.0;
		for (j=i+1;j<nspecies;j++)
		{
			pp[i+j*nspecies]= pp[j+i*nspecies] ; 
		}
	}
	parms->embedding_function = (double (*)(void *, double , double * ))exp_embedding; 
	parms->pass1_function = (EP (*)(void *, double )) exp_pass1; 
	parms->pass2_function = (EP (*)(void *, double )) exp_pass2; 
}
EP exp_pass1(EXP_PASS_PARMS*parms, double r2)
{
	EP ep;
	double r,ri; 
	ri = 1.0/sqrt(r2); 
	r = ri*r2; 

	ep.p = parms->f_e*exp(-parms->beta*(r*parms->r_e_inv - 1.0));
	ep.e = parms->phi_e*exp(-parms->gamma*(r*parms->r_e_inv - 1.0));
	return ep;
}
EP exp_pass2(EXP_PASS_PARMS*parms, double r2)
{
	EP dep;
	double r,ri; 

	ri = 1.0/sqrt(r2); 
	r = ri*r2; 
	dep.p = (-parms->beta*parms->r_e_inv*parms->f_e*exp(-parms->beta*(r*parms->r_e_inv - 1.0)))*ri;
	dep.e = (-parms->gamma*parms->r_e_inv*parms->phi_e*exp(-parms->gamma*(r*parms->r_e_inv - 1.0)))*ri;
	return dep;
}
double exp_embedding(EXP_EMBEDDING_PARMS*parms, double rho, double *dv_dr)
{
	double v,rr, x, y, T1, T2, T3, lnp, lnx,E_c,ab,gb,rho_e_inv;
	E_c = parms->E_c; 
	ab = parms->ab; 
	gb = parms->gb; 
	rho_e_inv = parms->rho_e_inv; 
	if (rho > 0.0) 
	{
		rr = rho*rho_e_inv;
/*
			x = pow(rr, ab);
			lnx = parms->ab*log(rr);
			y = pow(rr, gb);
			T1 = parms->E_c*x;
			T2 = parms->E_c*x*lnx;
			T3 = 6.0*parms->phi_e*y;
*/
		lnp = log(rr);
		y = exp(gb*lnp); 
		lnx = ab*lnp; 
		x = exp(lnx);

		T1 = x;
		T2 = x*lnx;
		T3 = y;
	
		v = E_c*(T2 -T1 - T3);
		*dv_dr = E_c*(ab*T2 - gb*T3)/rho;
	}
	else 
	{
		*dv_dr =0.0; 
		v=0.0; 
	}
	return v ;
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
