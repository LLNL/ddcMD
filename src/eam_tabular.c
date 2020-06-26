#include "eam_tabular.h"

#include "tfunction.h"
#include "ddcMalloc.h"
#include "species.h"
#include "mpiUtils.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))

typedef struct et_pass_parms_st
{
   TFUNCTION pairFunc;
} ET_PASS_PARMS;

typedef struct et_embedding_st
{
   TFUNCTION embedFunc;
} ET_EMBEDDING_PARMS;

static void et_embedMissingError(const char* species);
static void et_pairMissingError(const char* species1, const char* species2);

double et_embedding(ET_EMBEDDING_PARMS* parms, double rho, double* dv_dr)
{
   double v;
   tfunc_fdf(&parms->embedFunc, rho, &v, dv_dr);
   return v;
}

EP et_pass1(ET_PASS_PARMS* parms, double r2)
{
   EP ep;
   tfunc_f(&parms->pairFunc, r2, &ep.e, &ep.p);
   return ep;
}

EP et_pass2(ET_PASS_PARMS* parms, double r2)
{
   EP ep;
   tfunc_df(&parms->pairFunc, r2, &ep.e, &ep.p);
   return ep;
}

void eam_tabular_parms(POTENTIAL* object, EAM_PARMS* parms)
{
   parms->embedding_function = (double(*)(void*, double, double*)) et_embedding;
   parms->pass1_function = (EP (*)(void*, double))et_pass1; 
   parms->pass2_function = (EP (*)(void*, double))et_pass2;

   OBJECT* obj = (OBJECT*) object;
   
   unsigned nSpecies = parms->nspecies;    

   parms->pass_parms=ddcCalloc(nSpecies*nSpecies, sizeof(void*));
   ET_PASS_PARMS** pp = (ET_PASS_PARMS**) parms->pass_parms;
   for (unsigned ii=0; ii<nSpecies; ++ii) 
   {
      if (pp[ii+ii*nSpecies]==NULL)
	 ddcMallocAligned((void*)&pp[ii+ii*nSpecies],16,sizeof(ET_PASS_PARMS));
      for (unsigned jj=ii+1; jj<nSpecies; ++jj) 
      {
	 if (pp[ii+jj*nSpecies]==NULL)
	    ddcMallocAligned((void*)&pp[ii+jj*nSpecies],16,sizeof(ET_PASS_PARMS));
	 pp[jj+ii*nSpecies ] = pp[ii+jj*nSpecies];
      }
   }

   double rmax = 0.0;
   char** sName;
   species_get(NULL, NAMELIST, (void*)&sName);
   for (unsigned ii=0; ii<nSpecies; ++ii)
   {
      for (unsigned jj=ii; jj<nSpecies; ++jj)
      {
	 char keyword[256];
	 char* value;
	 snprintf(keyword, 255, "%s-%s_pair", sName[ii], sName[jj]);
	 if ( !object_testforkeyword(obj, keyword) )
	 {
	    snprintf(keyword, 255, "%s-%s_pair", sName[jj], sName[ii]);
	    if ( !object_testforkeyword(obj, keyword) )
	       et_pairMissingError(sName[ii], sName[jj]);
	 }
	 object_get(obj, keyword, &value, STRING, 1, "thisCantHappen");
	 pp[jj+ii*nSpecies]->pairFunc = tfunc_init(value);
//	 if (tfunc_rank(&pp[jj+ii*nSpecies]->pairFunc) != 2)
//	    et_wrongPairRankError(sName[ii], sName[jj]);
//	 rmax = MAX(rmax, tfunc_xMax(&pp[jj+ii*nSpecies]->pairFunc));
	 ddcFree(value);
      }
   }

   parms->embedding_parms = ddcCalloc(nSpecies, sizeof(ET_EMBEDDING_PARMS*));
   ET_EMBEDDING_PARMS** embedding_parms = (ET_EMBEDDING_PARMS**) parms->embedding_parms;
   for (unsigned ii=0; ii<nSpecies; ++ii)
   {
      embedding_parms[ii] = ddcCalloc(1, sizeof(ET_EMBEDDING_PARMS));
      char keyword[256];
      char* value;
      snprintf(keyword, 255, "%s_embed", sName[ii]);
      if ( !object_testforkeyword(obj, keyword))
	 et_embedMissingError(sName[ii]);
      object_get(obj, keyword, &value, STRING, 1, "thisCantHappen");
      embedding_parms[ii]->embedFunc = tfunc_init(value);
//      rmax = MAX(rmax, tfunc_xMax(&embedding_parms[ii]->embedFunc));
      ddcFree(value);
   }

   if (parms->rmax == 0.0)
      parms->rmax = rmax;

}


void et_embedMissingError(const char* species)
{
   if (getRank(0) == 0)
   {
      printf("ERROR:  Error in EAM tabular potential.\n"
	     "        No embedding function specified for species %s\n", species);
      abortAll(-1);
   }
}

void et_pairMissingError(const char* species1, const char* species2)
{
   if (getRank(0) == 0)
   {
      printf("ERROR:  Error in EAM tabular potential.\n"
	     "        No pair function specified for species pair %s-%s\n",
	     species1, species2);
      abortAll(-1);
   }
}

#ifdef COMPILE_UNUSED
static void et_wrongPairRankError(const char* species1, const char* species2)
{
   if (getRank(0) == 0)
   {
      printf("ERROR:  Error in EAM tabular potential.\n"
	     "        Pair function for species pair %s-%s does not have rank==2\n",
	     species1, species2);
      abortAll(-1);
   }
}
#endif



/* Local Variables: */
/* tab-width: 3 */
/* End: */
