#include "potential.h"
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "utilities.h"
#include "error.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"

void *mgpt_parms(POTENTIAL *);
void *eam_parms(POTENTIAL *);
void *eam1pass_parms(POTENTIAL *);
void *eam2pass_parms(POTENTIAL *);
void *pair_parms(POTENTIAL *);
void *charmm_parms(POTENTIAL *);
void *martini_parms(POTENTIAL *);
void *restraint_parms(POTENTIAL *);
void *eam_opt_parms(POTENTIAL *);
void *eam_onepass_parms(POTENTIAL *);
void *pairEnergy_parms(POTENTIAL *);
void *onebody_parms(POTENTIAL *);
void *reflect_parms(POTENTIAL *);
void *ewald_parms(POTENTIAL *);
void *fmm_parms(POTENTIAL *); //! FMM, Coulomb. Parameters
void *plasmaInit(POTENTIAL *);
void *orderSH_parms(POTENTIAL *);
//void *mirrorsym_parms(POTENTIAL*);
void *localYukawaInit(POTENTIAL *);
void *hycopInit(POTENTIAL *);
void *meam_parms(POTENTIAL*);
void *zeroPotential_parms(POTENTIAL *);

void mgpt(void *sys, void *parms, void *e);
void eam(void *sys, void *parms, void *e);
void eam1pass(void *sys, void *parms, void *e);
void eam2pass(void *sys, void *parms, void *e);
void pair(void *sys, void *parms, void *e);
void charmm(void *sys, void *parms, void *e);
void martini(void *sys, void *parms, void *e);
void restraint(void *sys, void *parms, void *e);
void eam_opt(void *sys, void *parms, void *e);
void eam_onepass(void *sys, void *parms, void *e);
void pairEnergy(void *sys, void *parms, void *e);
void onebody(void *sys, void *parms, void *e);
void reflect(void *sys, void *parms, void *e);
void ewald(void *sys, void *parms, void *e);
void plasma(void *sys, void *parms, void *e);
void orderSH(void *sys, void *parms, void *e);
//void mirrorsym(void *sys, void*parms, void *e);
void meam(void *sys, void*parms, void *e);
void zeroPotential(void *sys, void*parms, void *e);
void localYukawa(void *sys, void*parms, void *e);
void hycop(void *sys, void*parms, void *e);
void fmmCompute_ddcMD(void *sys, void *parms, void *e); // SYSTEM *system, FMM_PARMS *FMMParms, real *V

void plasma_write_dynamic(void *POTENTIAL, FILE *);

RCUT_TYPE *mgptCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *eamCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *eam1passCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *eam2passCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *pairCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *charmmCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *eam_optCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *eam_onepassCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *pairEnergyCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *onebodyCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *reflectCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *ewaldCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *plasmaCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *orderSHCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *fmmCutoff(void* sys, void *parms, int *n);
//RCUT_TYPE *mirrorsymCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *meamCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *zeroPotentialCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *restraintCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *localYukawaCutoff(void* sys, void *parms, int *n);
RCUT_TYPE *hycopCutoff(void* sys, void *parms, int *n);


POTENTIAL *potential_init(void *parent,char *name)
{
   POTENTIAL *potential;
   char *type;

   potential = (POTENTIAL *) object_initialize(name, "POTENTIAL", sizeof(POTENTIAL));
   potential->parent = parent;
   potential->itype = NO_POTENTIAL;
   potential->write_dynamics = NULL; 
   potential->getCutoffs = NULL;
   potential->eval_potential = NULL;
   potential->commMode = POT_TWOSIDED; // default potential assumes two sided ddc comms
   potential->use_gpu_list = 0;
   potential->neighborTableType = NEIGHBORTABLE_FAT; 
   object_get((OBJECT *) potential, "type", &type, STRING, 1, "NONE");
   potential->type = strdup(type); 
   if (strcmp(type, "MGPT") == 0)
   {
      potential->itype = MGPT;
      potential->parms = mgpt_parms(potential);
      potential->eval_potential = mgpt;
      potential->getCutoffs = mgptCutoff;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
   }
   if (strcmp(type, "REFLECT") == 0)
   {
      potential->itype = REFLECT;
      potential->parms = reflect_parms(potential);
      potential->eval_potential = reflect;
      potential->getCutoffs = reflectCutoff;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
   }
   if (strcmp(type, "ONEBODY") == 0)
   {
      potential->itype = ONEBODY;
      potential->parms = onebody_parms(potential);
      potential->eval_potential = onebody;
      potential->getCutoffs = onebodyCutoff;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
   }
   if (strcmp(type, "EAM") == 0)
   {
      potential->itype = EAM;
      potential->parms = eam_parms(potential);
      potential->eval_potential = eam;
      potential->getCutoffs = eamCutoff;
      potential->neighborTableType = NEIGHBORTABLE_SKINNY; 
   }
   if (strcmp(type, "EAM1PASS") == 0)
   {
      potential->itype = EAM1PASS;
      potential->parms = eam1pass_parms(potential);
      potential->eval_potential = eam1pass;
      potential->getCutoffs = eam1passCutoff;
      potential->neighborTableType = NEIGHBORTABLE_SKINNY; 
   }
   if (strcmp(type, "EAM2PASS") == 0)
   {
      potential->itype = EAM2PASS;
      potential->parms = eam2pass_parms(potential);
      potential->eval_potential = eam2pass;
      potential->getCutoffs = eam2passCutoff;
      potential->neighborTableType = NEIGHBORTABLE_SKINNY;   
   }
   if (strcmp(type, "EAM_OPT") == 0)
   {
      potential->itype = EAM_OPT;
      potential->parms = eam_opt_parms(potential);
      potential->eval_potential = eam_opt;
      potential->getCutoffs = eam_optCutoff;
      potential->neighborTableType = NEIGHBORTABLE_NONE; 
   }
   if (strcmp(type, "EAM_ONEPASS") == 0)
   {
      potential->itype = EAM_ONEPASS;
      potential->parms = eam_onepass_parms(potential);
      potential->eval_potential = eam_onepass;
      potential->getCutoffs = eam_onepassCutoff;
      potential->neighborTableType = NEIGHBORTABLE_NONE; 
   }
   if (strcmp(type, "PAIR") == 0)
   {
      potential->itype = PAIR;
      potential->eval_potential = pair;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
      potential->parms = pair_parms(potential);
      potential->getCutoffs = pairCutoff;
   }
   if (strcmp(type, "CHARMM") == 0)
   {
      potential->itype = CHARMM;
      potential->parms = charmm_parms(potential);
      potential->eval_potential = charmm;
      potential->getCutoffs = charmmCutoff;
      potential->neighborTableType = NEIGHBORTABLE_SKINNY; 
   }  
   if (strcmp(type, "MARTINI") == 0)
   {
      potential->itype = MARTINI;
      potential->eval_potential = martini;
      potential->neighborTableType = NEIGHBORTABLE_SKINNY; 
      potential->parms = martini_parms(potential);
      potential->getCutoffs = charmmCutoff;
   }  
   if (strcmp(type, "RESTRAINT") == 0)
   {
      potential->itype = RESTRAINT;
      potential->eval_potential = restraint;
      potential->getCutoffs = restraintCutoff;
      potential->neighborTableType = NEIGHBORTABLE_NONE; 
      potential->parms = restraint_parms(potential);
   }    
   if (strcmp(type, "LOCALYUKAWA") == 0)
   {
      potential->itype = LOCALYUKAWA;
      potential->parms = localYukawaInit(potential);
      potential->eval_potential = localYukawa;
      potential->getCutoffs = localYukawaCutoff;
      potential->neighborTableType = NEIGHBORTABLE_SKINNY; 
   }
   if (strcmp(type, "HYCOP") == 0)
   {
      potential->itype = HYCOP;
      potential->parms = hycopInit(potential);
      potential->eval_potential = hycop;
      potential->getCutoffs = hycopCutoff;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
   }
   if (strcmp(type, "PAIRENERGY") == 0)
   {
      potential->itype = PAIRENERGY;
      potential->parms = pairEnergy_parms(potential);
      potential->eval_potential = pairEnergy;
      potential->getCutoffs = pairEnergyCutoff;
      potential->neighborTableType = NEIGHBORTABLE_NONE; 
   }
   if (strcmp(type, "EWALD") == 0)
   {
      potential->itype = EWALD;
      potential->parms = ewald_parms(potential);
      potential->eval_potential = ewald;
      potential->getCutoffs = ewaldCutoff;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
   }
   if (strcmp(type, "PLASMA") == 0)
   {
      potential->itype = PLASMA;
      potential->parms = plasmaInit(potential);
      potential->write_dynamics = plasma_write_dynamic;  
      potential->getCutoffs = plasmaCutoff;
      potential->neighborTableType = NEIGHBORTABLE_NONE; 
   }
   if (strcmp(type, "ORDERSH") == 0)
   {
      potential->itype = ORDERSH;
      potential->parms = orderSH_parms(potential);
      potential->eval_potential = orderSH;
      potential->getCutoffs = orderSHCutoff;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
   }
   if (strcmp(type, "FMM") == 0) //! Fast multipole method. Coulomb 
   {
      potential->itype = FMM;
      potential->parms = fmm_parms(potential);
      potential->eval_potential = fmmCompute_ddcMD;
      potential->getCutoffs = fmmCutoff;
      potential->neighborTableType = NEIGHBORTABLE_NONE; 
      
   }
   if (strcmp(type, "NONE") == 0)
   {
      potential->itype = ZEROPOTENTIAL;
      potential->parms = zeroPotential_parms(potential);
      potential->eval_potential = zeroPotential;
      potential->getCutoffs = zeroPotentialCutoff;
      potential->neighborTableType = NEIGHBORTABLE_NONE; 
   }
/*
   if (strcmp(type, "MIRRORSYM") == 0)
   {
      potential->itype = MIRRORSYM;
      potential->parms = mirrorsym_parms(potential);
      potential->eval_potential = mirrorsym;
      potential->getCutoffs = mirrorsymCutoff;
      potential->neighborTableType = NEIGHBORTABLE_FAT; 
   }
*/
   if (strcmp(type, "MEAM") == 0)
   {
      potential->itype = MEAM;
      potential->parms = meam_parms(potential);
      potential->eval_potential = meam;
      potential->getCutoffs = meamCutoff;
      potential->neighborTableType = NEIGHBORTABLE_SKINNY; 
   }
   if (getRank(0) == 0 && potential->itype == NO_POTENTIAL) 
   {
      printf("ERROR:  Unspecified or unrecognized type in POTENTIAL object.\n"
        "        object name = %s\n"
        "        type = %s\n", potential->name, type);
      abortAll(-1);
   }
   assert(potential->eval_potential);
   assert(potential->getCutoffs);
   
   if (getRank(0) ==0) 
   {
      FILE* file = fopen("hpm.data","a"); 
      char* temp=strdup(potential->value); 
      if (strlen(temp) > 192 ) temp[193]='\0';
      fprintf(file,"******* %s ******\n",temp); 
      ddcFree(temp); 
      fclose(file);
   }
   return potential;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
