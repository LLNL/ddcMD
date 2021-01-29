#include "eam_at.h"

#include <math.h>
#include <mpi.h>

#include "ddcMalloc.h"
#include "utilities.h"
#include "species.h"
#include "units.h"
#include "mpiUtils.h"

/** Input parameters must be given in this order:
 *
 *  A alpha B b0 c c0 c1 c2 d
 *
 * i.e.,
 * Ta = 2.591061 1.05 91.0 2.8629 4.20 1.2157373 0.0271471 -0.1217350 4.076980;
 *
 * Units are:
 * A     (eV/Ang)
 * alpha (1/Ang)
 * B     (eV/Ang^3)
 * b0    (Ang)
 * c     (Ang)
 * c0    (eV/Ang^2)
 * c1    (eV/Ang^3)
 * c2    (eV/Ang^4)
 * d     (Ang)
 *
 */


typedef struct at_pass_parms_st
{
   double A;      // Embedding energy coefficient
   double alpha;  // Core correction length constant
   double B;      // Core correction coefficient
   double b0;     // Core correction cutoff
   double c;      // Pair function cutoff
   double c0;     // Pair function r^0 coefficient
   double c1;     // Pair function r^1 coefficient
   double c2;     // Pair function r^2 coefficient
   double d;      // Density function cutoff
   
} AT_PASS_PARMS;

typedef struct at_embedding_parms_st
{
   double negA;
} AT_EMBEDDING_PARMS;

double at_embedding(AT_EMBEDDING_PARMS* parms, double rho, double* dv_dr)
{
   double eps=1e-30; 
   double v; 
   v = parms->negA * sqrt(rho);
   *dv_dr = 0.5*v/(rho+eps);
   return v; 
}

EP at_pass1(AT_PASS_PARMS* parms, double r2)
{
   const double B = parms->B;
   const double b0 = parms->b0;
   const double alpha = parms->alpha;
   const double c = parms->c;
   const double c0 = parms->c0;
   const double c1 = parms->c1;
   const double c2 = parms->c2;
   const double d = parms->d;
   EP ep;
   double r = sqrt(r2);
   if (r>c)
      ep.e = 0;
   else
   {
      ep.e = (r-c)*(r-c) * (c0 + c1*r + c2*r2);
      if (r<b0)
	 ep.e += B*(b0-r)*(b0-r)*(b0-r)*exp(-alpha*r);
   }
   
   if (r>d)
      ep.p = 0;
   else
      ep.p = (r-d)*(r-d);

   return ep;
}

EP at_pass2(AT_PASS_PARMS* parms, double r2)
{
   const double B = parms->B;
   const double b0 = parms->b0;
   const double alpha = parms->alpha;
   const double c = parms->c;
   const double c0 = parms->c0;
   const double c1 = parms->c1;
   const double c2 = parms->c2;
   const double d = parms->d;
   
   EP ep; 
   double r = sqrt(r2);
   double ri = 1.0/r;
   if (r>c)
      ep.e = 0;
   else
   {
      ep.e = 2.0*(r-c)*(c0+c1*r+c2*r2) + (r-c)*(r-c)*(c1+2.0*c2*r);
      if (r<b0)
	 ep.e += -B*(b0-r)*(b0-r)*exp(-alpha*r) * (alpha*(b0-r)+3.0);
      ep.e *= ri;
   }

   if (r>d)
      ep.p = 0;
   else
      ep.p = 2.0*(r-d)*ri;
	
   return ep; 
}

void eam_at_parms(POTENTIAL* obj, EAM_PARMS* parms)
{
   
   parms->embedding_function = (double(*)(void*, double, double*)) at_embedding;
   parms->pass1_function = (EP (*)(void* , double))at_pass1; 
   parms->pass2_function = (EP (*)(void* , double))at_pass2;

   unsigned nSpecies = parms->nspecies;    
   if (nSpecies != 1 && getRank(0) == 0)
   {
      printf("Ackland-Thetford form supports single species only.\n"
	     "  Mixing rules are currently undefined.\n");
      MPI_Abort(MPI_COMM_WORLD,-1);
   }

   parms->pass_parms=ddcCalloc(nSpecies*nSpecies, sizeof(void*));
   AT_PASS_PARMS** pp = (AT_PASS_PARMS**) parms->pass_parms;
   for (unsigned ii=0; ii<nSpecies; ++ii) 
   {
      if (pp[ii+ii*nSpecies]==NULL)
	 ddcMallocAligned((void*)&pp[ii+ii*nSpecies],16,sizeof(AT_PASS_PARMS));
      for (unsigned jj=ii+1; jj<nSpecies; ++jj) 
      {
	 if (pp[ii+jj*nSpecies]==NULL)
	    ddcMallocAligned((void*)&pp[ii+jj*nSpecies],16,sizeof(AT_PASS_PARMS));
	 pp[jj+ii*nSpecies ] = pp[ii+jj*nSpecies];
      }
   }
   double length_convert = units_convert(1.0,"Angstrom",NULL); 
   double energy_convert = units_convert(1.0,"eV",NULL); 
   char** namelist;
   species_get(NULL,NAMELIST,(void *)&namelist);
   
   for (unsigned ii=0; ii<nSpecies; ++ii) 
   {
      AT_PASS_PARMS* pass_parms = pp[ii+ii*nSpecies];
      object_get((OBJECT *) obj, namelist[ii], &(pass_parms->A), DOUBLE, 9, "0.0"); 
      pass_parms->A     *= energy_convert/length_convert;
      pass_parms->alpha /= length_convert; 
      pass_parms->B     *= energy_convert/pow(length_convert, 3.0); 
      pass_parms->b0    *= length_convert; 
      pass_parms->c     *= length_convert; 
      pass_parms->c0    *= energy_convert/pow(length_convert, 2.0); 
      pass_parms->c1    *= energy_convert/pow(length_convert, 3.0); 
      pass_parms->c2    *= energy_convert/pow(length_convert, 4.0); 
      pass_parms->d     *= length_convert; 
   }

   parms->embedding_parms = ddcMalloc(nSpecies*sizeof(AT_EMBEDDING_PARMS*));
   for (unsigned ii=0; ii<nSpecies; ++ii)
   {
      parms->embedding_parms[ii] = ddcMalloc(sizeof(AT_EMBEDDING_PARMS));
      AT_EMBEDDING_PARMS* embedding_parms = parms->embedding_parms[ii];

      embedding_parms->negA = -pp[ii+ii*nSpecies]->A;
   }

   for (unsigned ii=0; ii<nSpecies-1; ++ii)
      for (unsigned jj=ii+1; jj<nSpecies; ++jj)
      {
	 //AT_PASS_PARMS* pass_parms = pp[jj+ii*nSpecies];
	 //define mixing rules here!
      }
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
