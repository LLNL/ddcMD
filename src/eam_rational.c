#include "eam_rational.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "object.h"

#include "species.h"
//#include "three_algebra.h"
#include "error.h"
#include "ptiming.h"
#include "ddcMalloc.h"
#include "units.h"

#include "mpiUtils.h"

static double rational_embedding(RATIONAL_EMBEDDING_PARMS *parms,double rho,double *dv_dr);
static void rational_pass0(RATIONAL_PASS_PARMS *parms,double r2,EP ep[1],EP dp[1]);
static EP   rational_pass1(RATIONAL_PASS_PARMS *parms,double r2);
static EP   rational_pass2(RATIONAL_PASS_PARMS *parms,double r2);

static ratexpr read_fit_object(OBJECT *fit) {
  double cutoff,*pcoeff,*qcoeff;
  int pdeg,qdeg;
  
  object_get(fit,"cutoff",&cutoff,DOUBLE,1,"0");
  object_get(fit,"orderP",&pdeg,INT,1,"0");
  object_get(fit,"orderQ",&qdeg,INT,1,"0");
  object_get(fit,"P",pcoeff = (double *) ddcCalloc(pdeg+1,sizeof(double)),DOUBLE,pdeg+1,"0.0");
  object_get(fit,"Q",qcoeff = (double *) ddcCalloc(qdeg+1,sizeof(double)),DOUBLE,qdeg+1,"0.0");

  /* Convert units of input parameters to simulation units */ {
    char *x_unit,*y_unit;
    object_get(fit,"xUnits",&x_unit,STRING,1,"NONE");
    object_get(fit,"yUnits",&y_unit,STRING,1,"NONE");
    const double x_conv = units_convert(1.0,x_unit,NULL);
    const double y_conv = units_convert(1.0,y_unit,NULL);    
    
    double f_p = y_conv,f_q = 1.0;

    for(int i = 0; i<=pdeg; i++) {
      pcoeff[i] *= f_p;
      f_p /= x_conv;
    }
    for(int i = 0; i<=qdeg; i++) {
      qcoeff[i] *= f_q;
      f_q /= x_conv;
    }

    cutoff *= x_conv;
  }

  /* Initialize array for streamlined/vectorizeable evaluation
     of function and derivative */
  double (*vd_p)[2] = (double (*)[2]) ddcCalloc(pdeg+1,sizeof(double [2]));
  double (*vd_q)[2] = (double (*)[2]) ddcCalloc(qdeg+1,sizeof(double [2]));
  
  vd_p[0][1] = 0.0;
  for(int i = 0; i<=pdeg; i++) {
    vd_p[i][0] = pcoeff[pdeg-i];
    if(i < pdeg)
      vd_p[i+1][1] = (pdeg-i)*pcoeff[pdeg-i];
  }
  
  vd_q[0][1] = 0.0;
  for(int i = 0; i<=qdeg; i++) {
    vd_q[i][0] = qcoeff[qdeg-i];
    if(i < qdeg)
      vd_q[i+1][1] = (qdeg-i)*qcoeff[qdeg-i];
  }

  if(getRank(0) == 0) {
    printf("\n\nP (deg = %d):\n",pdeg);
    for(int i = 0; i<=pdeg; i++)
      printf("i = %2d,  pcoeff = %10.3f   vd_p = %10.3f  %10.3f\n",
	     i,pcoeff[i],vd_p[i][0],vd_p[i][1]);
    printf("Q (deg = %d):\n",qdeg);
    for(int i = 0; i<=qdeg; i++)
      printf("i = %2d,  qcoeff = %10.3f   vd_q = %10.3f  %10.3f\n",
	     i,qcoeff[i],vd_q[i][0],vd_q[i][1]);
  }
  ratexpr rho = {
    cutoff ,
    { pdeg , pcoeff , vd_p } ,
    { qdeg , qcoeff , vd_q }
  };
  
  return rho;
}
//static double zero = 0.0,one = 1.0,zerozero[][2] = {{0.0,0.0}};
//static ratexpr zero_ratfun = { 0.0 , { 0 , &zero , zerozero } , { 0 , &one  , zerozero } };

void eam_rational_parms(POTENTIAL *object, EAM_PARMS *parms)
{
  const int nspecies = parms->nspecies;
  const char /* Suffixes used to generate names of: */
    F_fun_suffix[] = "_embedding", /* embedding (F), */
    rho_fun_suffix[] = "_density", /* density (rho), and */
    phi_fun_suffix[] = "_2body";   /* pair energy (phi) function objects */
  const char /* Different classes of density functions */
    rho_type_keyword[] = "density_type",
    rho_type_elementwise_keyword[]   = "elementwise"   , /* rho[i] = sum RHO_j( ||r_ij|| )  */
    rho_type_pairsymmetric_keyword[] = "pair_symmetric", /* RHO_ij(r) == RHO_ji(r)             */
    rho_type_pairgeneral_keyword[]   = "pair_general"  ; /* rho[i] = sum RHO_ji( ||r_ij || ) */


  RATIONAL_PASS_PARMS **pass_array = (RATIONAL_PASS_PARMS **)
    ddcCalloc(nspecies*nspecies,sizeof(RATIONAL_PASS_PARMS *));
  RATIONAL_EMBEDDING_PARMS **embed_array = (RATIONAL_EMBEDDING_PARMS **)
    ddcCalloc(nspecies,sizeof(RATIONAL_EMBEDDING_PARMS *));
  
  char **namelist;

  /* First half of pass_array contains debnsity (rho) functions,
     second half contains pair energy (phi) functions */

  for(int i = 0; i<nspecies; i++) {
    for(int j = 0; j<nspecies; j++) {
      pass_array[i*nspecies+j] = (RATIONAL_PASS_PARMS *) ddcCalloc(1,sizeof(RATIONAL_PASS_PARMS));
    }
    embed_array[i] = NULL;
  }

  parms->pass_parms = (void *) pass_array;
  parms->embedding_parms = (void *) embed_array;

  parms->embedding_function = (double (*)(void *,double,double *)) rational_embedding;
  parms->pass1_function = (EP (*)(void *,double)) rational_pass1;
  parms->pass2_function = (EP (*)(void *,double)) rational_pass2;

  species_get(NULL,NAMELIST,(void *) &namelist);
  /* Load embedding (F) functions */ {
    for (int i = 0; i<nspecies; i++) {
      OBJECT *fit;
      char fitname[80];
      
      snprintf(fitname,sizeof(fitname),"%s%s",namelist[i],F_fun_suffix);
      fit = object_initialize(fitname,"FIT",sizeof(OBJECT));
      
      if(fit == NULL) {
	fprintf(stderr,"Error @ %s:%d in %s(): Unable to initialize fit object %s.\n",
		__FILE__,__LINE__,__func__,fitname);
	exit(1);
      }
      
      embed_array[i] = (RATIONAL_EMBEDDING_PARMS *) ddcCalloc(1,sizeof(RATIONAL_EMBEDDING_PARMS));
      embed_array[i]->F_fun = read_fit_object(fit);    
    }
  }


  char *rho_type = NULL;
  object_get((OBJECT *) object,rho_type_keyword,&rho_type,STRING,1,"NONE");
  if(strcasecmp(rho_type,rho_type_elementwise_keyword) == 0) {

    for(int i = 0; i<nspecies; i++) {
      OBJECT *fit;
      char fitname[80];
      
      snprintf(fitname,sizeof(fitname),"%s%s",namelist[i],rho_fun_suffix);
      fit = object_initialize(fitname,"FIT",sizeof(OBJECT));
      
      if(fit == NULL) {
	fprintf(stderr,"Error @ %s:%d in %s(): Unable to initialize fit object %s.\n",
		__FILE__,__LINE__,__func__,fitname);
	exit(1);
      }
    
      {
	ratexpr fun = read_fit_object(fit);
	for(int j = 0; j<nspecies; j++)
	  pass_array[j*nspecies + i]->rho_fun = fun;
      }
    }

  } else if(strcasecmp(rho_type,rho_type_pairsymmetric_keyword) == 0) {

    for(int i = 0; i<nspecies; i++)
      for(int j = i; j<nspecies; j++) {
	OBJECT *fit;
	char fitname[80];
      
	/* Initialize density function object */ {
	  snprintf(fitname,sizeof(fitname),"%s_%s%s",namelist[i],namelist[j],rho_fun_suffix);
	  fit = object_initialize(fitname,"FIT",sizeof(OBJECT));
	  if(fit == NULL) {
	    char fitname2[80];
	    snprintf(fitname,sizeof(fitname2),"%s_%s%s",namelist[j],namelist[i],rho_fun_suffix);
	    fit = object_initialize(fitname2,"FIT",sizeof(OBJECT));
	    if(fit == NULL) {
	      fprintf(stderr,"Error @ %s:%d in %s(): Unable to initialize fit object %s or %s.\n",
		      __FILE__,__LINE__,__func__,fitname,fitname2);
	      exit(1);
	    }
	  }
	}

	/* Allocate entry and insert into pair array */ {
	  ratexpr fun = read_fit_object(fit);
	  pass_array[i*nspecies+j]->rho_fun = fun;
	  pass_array[j*nspecies+i]->rho_fun = fun;
	}
      }
    
    } else if(strcasecmp(rho_type,rho_type_pairgeneral_keyword) == 0) {
      
      for(int i = 0; i<nspecies; i++)
	for(int j = 0; j<nspecies; j++) {
	  OBJECT *fit;
	  char fitname[80];
	  
	  snprintf(fitname,sizeof(fitname),"%s_%s%s",namelist[i],namelist[j],rho_fun_suffix);
	  fit = object_initialize(fitname,"FIT",sizeof(OBJECT));
	  if(fit == NULL) {
	    fprintf(stderr,"Error @ %s:%d in %s(): Unable to initialize fit object %s.\n",
		    __FILE__,__LINE__,__func__,fitname);
	    exit(1);
	  }
	  
	  /* Allocate entry and insert into pair array */ {
	    ratexpr fun = read_fit_object(fit);
	    pass_array[i*nspecies+j]->rho_fun = fun;
	  }
	}
      
    } else {
      fprintf(stderr,
	      "Error @ %s:%d in %s(): Unable to understand rho_type='%s', "
	      "should be 'elementwise', 'pairsymmetric', or 'pairgeneral'.\n",
	      __FILE__,__LINE__,__func__,rho_type);
      exit(1);
    }


    /* Load pair energy functions */ {

      for(int i = 0; i<nspecies; i++)
	for(int j = i; j<nspecies; j++) {
	  OBJECT *fit;
	  char fitname[80];
	  
	  /* Initialize density function object */ {
	    snprintf(fitname,sizeof(fitname),"%s_%s%s",namelist[i],namelist[j],phi_fun_suffix);
	    fit = object_initialize(fitname,"FIT",sizeof(OBJECT));
	    if(fit == NULL) {
	      char fitname2[80];
	      snprintf(fitname2,sizeof(fitname2),"%s_%s%s",namelist[j],namelist[i],phi_fun_suffix);
	      fit = object_initialize(fitname2,"FIT",sizeof(OBJECT));
	      if(fit == NULL) {
		fprintf(stderr,"Error @ %s:%d in %s(): Unable to initialize fit object %s or %s.\n",
			__FILE__,__LINE__,__func__,fitname,fitname2);
		exit(1);
	      }
	    }
	  }
	  
	  /* Allocate entry and insert into pair array */ {
	    ratexpr fun = read_fit_object(fit);
	    pass_array[i*nspecies+j]->phi_fun = fun;
	    pass_array[j*nspecies+i]->phi_fun = fun;
	  }
	}
    }
}

static double eval_poly(int deg,const double coeff[],double x) {
  double y = 0.0;
  if(deg >= 0) {
    int i;
    y = coeff[deg];
    for(i = deg-1; i>= 0; i--)
      y = coeff[i] + x*y;
  }
  return y;
}
static double eval_poly_deriv(int deg,const double coeff[],double x) {
  double dy = 0.0;
  if(deg >= 1) {
    int i;
    double ir = deg-1;;
    dy = deg*coeff[deg];
    for(i = deg-1; i>= 1; i--) {
      dy = ir*coeff[i] + x*dy;
      ir -= 1.0;
    }
  }
  return dy;
}
double eval_rational(int num_deg,const double num_coeff[],
		     int den_deg,const double den_coeff[],
		     double x,double deriv[1]) {
  const double p = eval_poly(num_deg,num_coeff,x);
  const double qinv = 1.0/eval_poly(den_deg,den_coeff,x);
  const double fun_val = p*qinv;

#if 0
  printf("eval rat. pdeg=%d qdeg=%d p=0x%llx q=0x%8llx  "
	 "x = %20.10f  fun_val=%20.10f\n",
	 num_deg,den_deg,
	 (unsigned long long int) num_coeff,
	 (unsigned long long int) den_coeff,
	 x,
	 fun_val);
#endif

  if(deriv != NULL) {
    const double dp = eval_poly_deriv(num_deg,num_coeff,x);
    const double dq = eval_poly_deriv(den_deg,den_coeff,x);
    deriv[0] = qinv*(dp - fun_val*dq);
  }
  return fun_val;
}


double rational_embedding(RATIONAL_EMBEDDING_PARMS *parms, double rho, double dF_drho[1])
{
  double F = 0.0;

  /* Evaluate embedding function */
  if(rho < parms->F_fun.cutoff) {
    const int pdeg = parms->F_fun.numerator.degree;
    const int qdeg = parms->F_fun.denominator.degree;
    const double *pp = parms->F_fun.numerator.coeff;
    const double *qq = parms->F_fun.denominator.coeff;
    
    /* if dF_drho == NULL, then derivative does not get evaluated. */
    F = eval_rational(pdeg,pp,qdeg,qq,rho,dF_drho);
  } else if(dF_drho != NULL)
    dF_drho[0] = 0.0;

  return F;
}

void rational_pass0(RATIONAL_PASS_PARMS *parms,double r2,EP *ep,EP *dp)
{
  //const double ri = 1.0/sqrt(r2),r = r2*ri;

  /* Evaluate rho function */
  if(r2 < parms->rho_fun.cutoff) {
    const int pdeg = parms->rho_fun.numerator.degree;
    const int qdeg = parms->rho_fun.denominator.degree;
    const double *pp = parms->rho_fun.numerator.coeff;
    const double *qq = parms->rho_fun.denominator.coeff;
    
    if(dp != NULL)
      ep->p = eval_rational(pdeg,pp,qdeg,qq,r2,&(dp->p));
    else
      ep->p = eval_rational(pdeg,pp,qdeg,qq,r2,NULL);
  } else {
    ep->p = 0.0;
    if(dp != NULL) dp->p = 0.0;
  }
  
  /* Evaluate pair function */
  if(r2 < parms->phi_fun.cutoff) {
    const int pdeg = parms->phi_fun.numerator.degree;
    const int qdeg = parms->phi_fun.denominator.degree;
    const double *pp = parms->phi_fun.numerator.coeff;
    const double *qq = parms->phi_fun.denominator.coeff;
    
    if(dp != NULL)
      ep->e = eval_rational(pdeg,pp,qdeg,qq,r2,&(dp->e));
    else
      ep->e = eval_rational(pdeg,pp,qdeg,qq,r2,NULL);
  } else {
    ep->e = 0.0;
    if(dp != NULL) dp->e = 0.0;
  }
  
  /* What we want to return is -(1/r) d{rho,pair}(r^2)/dr */
  if(dp != NULL) {
    dp->p *= 2.0;
    dp->e *= 2.0;
  }
  
}

EP rational_pass1(RATIONAL_PASS_PARMS *parms, double r2)
{
  EP ep; 
  rational_pass0(parms,r2,&ep,NULL);
  return ep; 
}

EP rational_pass2(RATIONAL_PASS_PARMS *parms, double r2)
{
  EP ep,dp; 
  rational_pass0(parms,r2,&ep,&dp);
  return dp; 
}
