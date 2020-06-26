#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include "analysis.h"
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "mpiUtils.h"
#include "simulate.h"

#include "npartHist.h"

void *kineticEnergyDistn_parms(ANALYSIS *);
void kineticEnergyDistn_output(ANALYSIS *analysis);
void kineticEnergyDistn_eval(ANALYSIS *analysis);
void kineticEnergyDistn_clear(ANALYSIS *analysis);

void *coarsegrain_parms(ANALYSIS *);
void coarsegrain_output(ANALYSIS *analysis);
void coarsegrain_eval(ANALYSIS *analysis);
void coarsegrain_clear(ANALYSIS *analysis);

void *subsetWrite_parms(ANALYSIS *);
void subsetWrite(ANALYSIS *analysis);
void subsetWrite_eval(ANALYSIS *analysis);

void *velocityAutocorrelation_parms(ANALYSIS *analysis);
void velocityAutocorrelation_output(ANALYSIS *analysis);
void velocityAutocorrelation_eval(ANALYSIS *analysis);

void *plasmaWake_parms(ANALYSIS *analysis);
void plasmaWake_output(ANALYSIS *analysis);
void plasmaWake_eval(ANALYSIS *analysis);
void plasmaWake_clear(ANALYSIS *analysis);

void *forceAverage_parms(ANALYSIS *analysis);
void forceAverage_output(ANALYSIS *analysis);
void forceAverage_eval(ANALYSIS *analysis);
void forceAverage_clear(ANALYSIS *analysis);

void *paircorrelation_parms(ANALYSIS *analysis);
void paircorrelation_output(ANALYSIS *analysis);
void paircorrelation_eval(ANALYSIS *analysis);
void paircorrelation_clear(ANALYSIS* analysis);

void *pairAnalysis_parms(ANALYSIS *analysis);
void pairAnalysis_output(ANALYSIS *analysis);
void pairAnalysis_eval(ANALYSIS *analysis);
void pairAnalysis_clear(ANALYSIS* analysis);

void *stressWrite_parms(ANALYSIS *);
void stressWrite_output(ANALYSIS *analysis);
void stressWrite_eval(ANALYSIS *analysis);

void *ssf_parms(ANALYSIS *);
void ssf_output(ANALYSIS *analysis);
void ssf_eval(ANALYSIS *analysis);

void *dsf_parms(ANALYSIS *);
void dsf_output(ANALYSIS *analysis);
void dsf_eval(ANALYSIS *analysis);

void *vcmWrite_parms(ANALYSIS *);
void vcmWrite_output(ANALYSIS *analysis);
void vcmWrite_eval(ANALYSIS *analysis);

void* centroSym_parms(ANALYSIS*);
void centroSym_output(ANALYSIS* analysis);
void centroSym_eval(ANALYSIS* analysis);

void* AcklandJones_parms(ANALYSIS*);
void AcklandJones_output(ANALYSIS* analysis);
void AcklandJones_eval(ANALYSIS* analysis);

void* vonMises_parms(ANALYSIS*);
void vonMises_output(ANALYSIS* analysis);
void vonMises_eval(ANALYSIS* analysis);

void* quaternion_parms(ANALYSIS*);
void quaternion_output(ANALYSIS* analysis);
void quaternion_eval(ANALYSIS* analysis);

void* boundState_parms(ANALYSIS*);
void boundState_output(ANALYSIS* analysis);
void boundState_eval(ANALYSIS* analysis);

void* nearMiss_parms(ANALYSIS*);
void nearMiss_output(ANALYSIS* analysis);
void nearMiss_eval(ANALYSIS* analysis);

void* zdensity_parms(ANALYSIS*);
void zdensity_output(ANALYSIS* analysis);
void zdensity_eval(ANALYSIS* analysis);

void* ewaldTransportCoeff_parms(ANALYSIS*);
void ewaldTransportCoeff_output(ANALYSIS* analysis);
void ewaldTransportCoeff_eval(ANALYSIS* analysis);

void*localYukawaAnalysis_parms(ANALYSIS*);
void localYukawaAnalysis_output(ANALYSIS* analysis);
void localYukawaAnalysis_eval(ANALYSIS* analysis);

void* ewaldAnalysis_parms(ANALYSIS*);
void  ewaldAnalysis_output(ANALYSIS* analysis);
void  ewaldAnalysis_eval(ANALYSIS* analysis);

void* dataSubset_parms(ANALYSIS*);
void  dataSubset_output(ANALYSIS* analysis);
void  dataSubset_eval(ANALYSIS* analysis);
void  dataSubset_startup(ANALYSIS* analysis);

void* projectileRV_parms(ANALYSIS*);
void  projectileRV_output(ANALYSIS* analysis);
void  projectileRV_eval(ANALYSIS* analysis);
void  projectileRV_clear(ANALYSIS* analysis);

void* event_parms(ANALYSIS*);
void  event_eval(ANALYSIS*);
void  event_output(ANALYSIS*);

void* cdb_parms(ANALYSIS*);
void  cdb_eval(ANALYSIS*);
void  cdb_output(ANALYSIS*);

void* cholAnalysis_parms(ANALYSIS*);
void  cholAnalysis_eval(ANALYSIS*);
void  cholAnalysis_output(ANALYSIS*);

void analysis_NULL(ANALYSIS* analysis)
{
}
void analysis_startup(ANALYSIS* analysis)
{
 
 SIMULATE *simulate=analysis->parent; 
 if ( TEST0(simulate->loop, analysis->eval_rate) ) analysis->eval(analysis);
  analysis->clear(analysis);
}

ANALYSIS *analysis_init(void *parent,char *name)
{
	ANALYSIS *analysis;
	char *type;

	analysis = (ANALYSIS *) object_initialize(name, "ANALYSIS", sizeof(ANALYSIS));
	analysis->parent = parent;
	analysis->itype = NO_ANALYSIS;
	analysis->clear = analysis_NULL;
	analysis->startup = analysis_startup;
	object_get((OBJECT *) analysis, "type", &type, STRING, 1, "NONE");
	object_get((OBJECT *) analysis, "eval_rate", &analysis->eval_rate, INT, 1, "0");
	object_get((OBJECT *) analysis, "outputrate", &analysis->outputrate, INT, 1, "0");
	if (strncasecmp(type, "KINETICENERGYDISTN", strlen("KINETICENERGYDISTN")) == 0)
	{
	   analysis->itype = KINETICENERGYDISTN;
	   analysis->parms = kineticEnergyDistn_parms(analysis);
	   analysis->output= kineticEnergyDistn_output;
	   analysis->eval  = kineticEnergyDistn_eval;
	   analysis->clear = kineticEnergyDistn_clear;
	}
	if (strncasecmp(type, "COARSEGRAIN", strlen("COARSEGRAIN")) == 0)
	{
	   analysis->itype = COARSEGRAIN;
	   analysis->parms = coarsegrain_parms(analysis);
	   analysis->output = coarsegrain_output;
	   analysis->eval = coarsegrain_eval;
	   analysis->clear = coarsegrain_clear;
	}
	if (strcasecmp(type, "subset_write") == 0 ||
	    strcasecmp(type, "subsetWrite") == 0 )
	{
	   analysis->itype = SUBSET_WRITE;
	   analysis->parms = subsetWrite_parms(analysis);
	   analysis->output = subsetWrite;
	   analysis->eval = subsetWrite_eval;
	}
	if (strncasecmp(type, "VELOCITYAUTOCORRELATION", strlen("VELOCITYAUTOCORRELATION")) == 0)
	{
	   analysis->itype = VELOCITYAUTOCORRELATION;;
	   analysis->parms = velocityAutocorrelation_parms(analysis);
	   analysis->output = velocityAutocorrelation_output;
	   analysis->eval = velocityAutocorrelation_eval;
	}
	if (strncasecmp(type, "PAIRCORRELATION", strlen("PAIRCORRELATION")) == 0)
	{
	   analysis->itype = PAIRCORRELATION;
	   analysis->parms = paircorrelation_parms(analysis);
	   analysis->output = paircorrelation_output;
	   analysis->eval = paircorrelation_eval;
	   analysis->clear = paircorrelation_clear;
	}
	if (strncasecmp(type, "PAIRANALYSIS", strlen("PAIRANALYSIS")) == 0)
	{
	   analysis->itype = PAIRANALYSIS;
	   analysis->parms = pairAnalysis_parms(analysis);
	   analysis->output = pairAnalysis_output;
	   analysis->eval = pairAnalysis_eval;
	   analysis->clear = pairAnalysis_clear;
	}
	if (strncasecmp(type, "PLASMA_WAKE", strlen("PLASMA_WAKE")) == 0)
	{
	   analysis->itype = PLASMAWAKE;
	   analysis->parms = plasmaWake_parms(analysis);
	   analysis->output = plasmaWake_output;
	   analysis->eval = plasmaWake_eval;
	   analysis->clear = plasmaWake_clear;
	}
	if (strncasecmp(type, "FORCE_AVERAGE", strlen("FORCE_AVERAGE")) == 0)
	{
	   analysis->itype = FORCEAVERAGE;
	   analysis->parms = forceAverage_parms(analysis);
	   analysis->output = forceAverage_output;
	   analysis->eval = forceAverage_eval;
	   analysis->clear = forceAverage_clear;
	}
	if (strncasecmp(type, "STRESS_WRITE", strlen("STRESS_WRITE")) == 0)
	{
	   analysis->itype = STRESS_WRITE;
	   analysis->parms = stressWrite_parms(analysis);
	   analysis->output = stressWrite_output;
	   analysis->eval = stressWrite_eval;
	}
	if (strncasecmp(type, "SSF", strlen("SSF")) == 0 )
	{
	   analysis->itype = STATIC_STRUCTURE_FACTOR;
	   analysis->parms = ssf_parms(analysis);
	   analysis->output = ssf_output;
	   analysis->eval = ssf_eval;
	}
	if (strncasecmp(type, "DSF", strlen("DSF")) == 0 ||
	    strncasecmp(type, "DynamicStructureFactor", strlen("DynamicStructureFactor")) == 0 ||
	    strncasecmp(type, "Dynamic_Structure_Factor", strlen("Dynamic_Structure_Factor")) == 0 )
	{
	   analysis->itype = DYNAMIC_STRUCTURE_FACTOR;
	   analysis->parms = dsf_parms(analysis);
	   analysis->output = dsf_output;
	   analysis->eval = dsf_eval;
	}
	if (strcasecmp(type, "vcmWrite") == 0 ||
	    strcasecmp(type, "vcm_write") == 0 )
	{
	   analysis->itype = VCM_WRITE;
	   analysis->parms = vcmWrite_parms(analysis);
	   analysis->output = vcmWrite_output;
	   analysis->eval = vcmWrite_eval;
	}

	if (strcasecmp(type, "vcmWrite") == 0 ||
	    strcasecmp(type, "vcm_write") == 0 )
	{
	   analysis->itype = VCM_WRITE;
	   analysis->parms = vcmWrite_parms(analysis);
	   analysis->output = vcmWrite_output;
	   analysis->eval = vcmWrite_eval;
	}
	if (strcasecmp(type, "centroSymmetry") == 0 ||
	    strcasecmp(type, "CENTRO_SYMMETRY") == 0 ||
	    strcasecmp(type, "centralSymmetry") == 0 ||
	    strcasecmp(type, "CENTRAL_SYMMETRY") == 0 )
	{
	   analysis->itype = CENTRO_SYMMETRY;
	   analysis->parms = centroSym_parms(analysis);
	   analysis->output = centroSym_output;
	   analysis->eval = centroSym_eval;
	}
	if (strcasecmp(type, "acklandJones") == 0 ||
	    strcasecmp(type, "ACKLAND_JONES") == 0 )
	{
	   analysis->itype = ACKLAND_JONES;
	   analysis->parms = AcklandJones_parms(analysis);
	   analysis->output = AcklandJones_output;
	   analysis->eval = AcklandJones_eval;
	}
	if (strcasecmp(type, "vonMises") == 0 ||
	    strcasecmp(type, "VON_MISES") == 0 )
	{
	   analysis->itype = VON_MISES;
	   analysis->parms = vonMises_parms(analysis);
	   analysis->output = vonMises_output;
	   analysis->eval = vonMises_eval;
	}
	if (strcasecmp(type, "quaternion") == 0 ||
	    strcasecmp(type, "QUATERNION") == 0 )
	{
	   analysis->itype = QUATERNION;
	   analysis->parms = quaternion_parms(analysis);
	   analysis->output = quaternion_output;
	   analysis->eval = quaternion_eval;
	}
	if (strcasecmp(type, "boundState") == 0 ||
	    strcasecmp(type, "BOUND_STATE") == 0 )
	{
	   analysis->itype = BOUND_STATE;
	   analysis->parms = boundState_parms(analysis);
	   analysis->output = boundState_output;
	   analysis->eval = boundState_eval;
	}

	if (strcasecmp(type, "nearMiss") == 0 ||
	    strcasecmp(type, "NEAR_MISS") == 0 )
	{
	   analysis->itype = NEAR_MISS;
	   analysis->parms = nearMiss_parms(analysis);
	   analysis->output = nearMiss_output;
	   analysis->eval = nearMiss_eval;
	}

	if (strcasecmp(type, "ewaldTransportCoeff") == 0  ||
	    strcasecmp(type, "EWALD_TRANSPORT_COEFF") == 0 )
	{
	   analysis->itype = EWALD_TRANSPORT_COEFF;
	   analysis->parms = ewaldTransportCoeff_parms(analysis);
	   analysis->output = ewaldTransportCoeff_output;
	   analysis->eval = ewaldTransportCoeff_eval;
	}
	if (strcasecmp(type, "zdensity") == 0 )
	{
	   analysis->itype = ZDENSITY;
	   analysis->parms = zdensity_parms(analysis);
	   analysis->output = zdensity_output;
	   analysis->eval = zdensity_eval;
	}
	if (strcasecmp(type, "DATASUBSET") == 0 )
	{
	   analysis->itype = DATASUBSET;
	   analysis->parms = dataSubset_parms(analysis);
	   analysis->output = dataSubset_output;
	   analysis->eval = dataSubset_eval;
	   analysis->startup = dataSubset_startup;
	}

	if (strcasecmp(type, "PROJECTILE_RV") == 0 )
	{
	   analysis->itype = PROJECTILE_RV;
	   analysis->parms = projectileRV_parms(analysis);
	   analysis->output = projectileRV_output;
	   analysis->eval = projectileRV_eval;
      analysis->clear = projectileRV_clear;
	}

	if(strcasecmp(type,"particle_histogram") == 0) 
   {
	  analysis->itype  = PARTICLE_HISTOGRAM;
	  analysis->parms  = npartHist_parms(analysis);
	  analysis->output = npartHist_output;
	  analysis->eval   = npartHist_eval;
	  analysis->clear  = npartHist_clear;
	  //analysis->startup = npartHist_clear;
	}

	if(strcasecmp(type,"localYukawaAnalysis") == 0) 
   {
	  analysis->itype = LOCALYUKAWA_ANALYSIS;
	  analysis->parms  = localYukawaAnalysis_parms(analysis);
	  analysis->output = localYukawaAnalysis_output;
	  analysis->eval   = localYukawaAnalysis_eval;
	  analysis->startup= analysis_NULL;
	}
	if(strcasecmp(type,"ewaldAnalysis") == 0) 
   {
	  analysis->itype = EWALD_ANALYSIS;
	  analysis->parms  = ewaldAnalysis_parms(analysis);
	  analysis->output = ewaldAnalysis_output;
	  analysis->eval   = ewaldAnalysis_eval;
	}
	if(strcasecmp(type,"events") == 0) 
   {
	  analysis->itype = EVENTS;
	  analysis->parms  = event_parms(analysis);
	  analysis->output = event_output;
  	  analysis->eval   = event_eval;
 }
   if(strcasecmp(type,"cdb") == 0) 
   {
      analysis->itype = CDB_ANALYSIS;
      analysis->parms  = cdb_parms(analysis);
      analysis->output = cdb_output;
      analysis->eval   = cdb_eval;
      analysis->startup  = cdb_eval;
   }
   if(strcasecmp(type,"cholAnalysis") == 0) 
   {
      analysis->itype = CHOL_ANALYSIS;
      analysis->parms  = cholAnalysis_parms(analysis);
      analysis->output = cholAnalysis_output;
      analysis->eval   = cholAnalysis_eval;
      analysis->startup  = cholAnalysis_eval;
   }

   if (analysis->itype == NO_ANALYSIS && getRank(0) == 0)
   {
      printf("Error: Unsupported analysis type %s\n Game Over\n", type);
      exit(1);
   }
   return analysis;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
