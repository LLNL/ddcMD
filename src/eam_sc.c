#include "eam_sc.h"

#include <math.h>
#include <mpi.h>

#include "ddcMalloc.h"
#include "utilities.h"
#include "species.h"
#include "units.h"
#include "mpiUtils.h"

/** Input parameters must be given in this order:
 *
 *  a(angstroms), c, epsilon(eV), n, m
 *
 * i.e., Fe = 3.4714 24.939 0.0173 8.137 4.788;
 *
 */


typedef struct sc_pass_parms_st
{
   double a;      // length units (input Angstroms)
   double c;
   double epsilon; // energy units (input eV)
   double n;
   double m;
   double a2;
   double half_n;
   double half_m;
} SC_PASS_PARMS;

typedef struct sc_embedding_parms_st
{
   double neg_c_times_epsilon;
} SC_EMBEDDING_PARMS;

double sc_embedding(SC_EMBEDDING_PARMS* parms, double rho, double* dv_dr)
{
   double v; 
   v = parms->neg_c_times_epsilon * sqrt(rho);
   *dv_dr = 0.5*v/rho;
   return v; 
}

void sc_pass0(SC_PASS_PARMS* parms, double r2, EP* ep,  EP* dp)
{
   double arg2 = parms->a2/r2;
   double ri  = 1.0/sqrt(r2);
   ep->e = parms->epsilon * pow(arg2, parms->half_n);
   ep->p = pow(arg2, parms->half_m);
   dp->e = -parms->n * ep->e * ri;
   dp->p = -parms->m * ep->p * ri;
}

EP sc_pass1(SC_PASS_PARMS* parms, double r2)
{
   EP ep;
/*    double arg2 = parms->a2/r2; */
/*    ep.e = parms->epsilon * pow(arg2, parms->half_n); */
/*    ep.p =                  pow(arg2, parms->half_m); */

   double ri = 1.0/sqrt(r2);
   double arg = parms->a*ri;
   ep.e = parms->epsilon * pow(arg, parms->n);
   ep.p =                  pow(arg, parms->m);


   return ep;
}

EP sc_pass2(SC_PASS_PARMS* parms, double r2)
{
   EP ep; 
   double arg2 = parms->a2/r2;
   ep.e = -(parms->n * parms->epsilon * pow(arg2, parms->half_n))/r2;
   ep.p = -(parms->m *                  pow(arg2, parms->half_m))/r2;
   return ep; 
}

void eam_sc_parms(POTENTIAL* obj, EAM_PARMS* parms)
{
   
   parms->embedding_function = (double(*)(void*, double, double*)) sc_embedding;
   parms->pass1_function = (EP (*)(void* , double))sc_pass1; 
   parms->pass2_function = (EP (*)(void* , double))sc_pass2;

   unsigned nSpecies = parms->nspecies;    
   if (nSpecies != 1 && getRank(0) == 0)
   {
      printf("Sutton-Chen form supports single species only.\n"
	     "  Mixing rules are currently undefined.\n");
      MPI_Abort( MPI_COMM_WORLD,-1);
   }

   parms->pass_parms=ddcCalloc(nSpecies*nSpecies, sizeof(void*));
   SC_PASS_PARMS** pp = (SC_PASS_PARMS**) parms->pass_parms;
   for (unsigned ii=0; ii<nSpecies; ++ii) 
   {
      if (pp[ii+ii*nSpecies]==NULL)
	 ddcMallocAligned((void*)&pp[ii+ii*nSpecies],16,sizeof(SC_PASS_PARMS));
      for (unsigned jj=ii+1; jj<nSpecies; ++jj) 
      {
	 if (pp[ii+jj*nSpecies]==NULL)
	    ddcMallocAligned((void*)&pp[ii+jj*nSpecies],16,sizeof(SC_PASS_PARMS));
	 pp[jj+ii*nSpecies ] = pp[ii+jj*nSpecies];
      }
   }
   double length_convert = units_convert(1.0,"Angstrom",NULL); 
   double energy_convert = units_convert(1.0,"eV",NULL); 
   char** namelist;
   species_get(NULL,NAMELIST,(void *)&namelist);
   
   for (unsigned ii=0; ii<nSpecies; ++ii) 
   {
      SC_PASS_PARMS* pass_parms = pp[ii+ii*nSpecies];
      object_get((OBJECT *) obj, namelist[ii], &(pass_parms->a), DOUBLE, 5, "0.0"); 
      pass_parms->a *= length_convert; 
      pass_parms->epsilon *= energy_convert;
      pass_parms->a2 = pass_parms->a * pass_parms->a;
      pass_parms->half_m = 0.5*pass_parms->m;
      pass_parms->half_n = 0.5*pass_parms->n;
   }

   parms->embedding_parms = ddcMalloc(nSpecies*sizeof(SC_EMBEDDING_PARMS*));
   for (unsigned ii=0; ii<nSpecies; ++ii)
   {
      double epsilon = pp[ii+ii*nSpecies]->epsilon;
      double c = pp[ii+ii*nSpecies]->c;
      parms->embedding_parms[ii] = ddcMalloc(sizeof(SC_EMBEDDING_PARMS));
      SC_EMBEDDING_PARMS* embedding_parms = parms->embedding_parms[ii];
      embedding_parms->neg_c_times_epsilon = -c*epsilon;
   }

   for (unsigned ii=0; ii<nSpecies-1; ++ii)
      for (unsigned jj=ii+1; jj<nSpecies; ++jj)
      {
	 SC_PASS_PARMS* pass_parms = pp[jj+ii*nSpecies];
	 //define mixing rules here!
	 pass_parms->a2 = pass_parms->a * pass_parms->a;
	 pass_parms->half_m = 0.5*pass_parms->m;
	 pass_parms->half_n = 0.5*pass_parms->n;
      }
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
