#ifndef HYCOP_H
#define HYCOP_H
#include "eq.h"
#include "potential.h"
#include "group.h"
#include "random.h"
#include "system.h"
#include "species.h"
#include "langevin.h"
#include "expandbuffer.h"
typedef struct hycopLangevin_parms_st
{
   LANGEVIN_PARMS_COMMON

	EQTARGET *Teq_vs_time;
	double Teq; 

   POTENTIAL *potential; 
   int nspecies; 
   SPECIES **species; 
   double **rate;
   double *D; 
} HYCOPLANGEVIN_PARMS;
typedef struct hycopPotential_parms_st
{
   
	POTENTIAL *potential; 
   SYSTEM *system; 
	unsigned nspecies; 
} HYCOP_PARMS ; 
void hycopLangevin_Update1(GROUP *group, int mode, STATE *state, double time, double dt);
void hycopLangevin_start(GROUP *group, int mode, STATE *state, double time, double dt);
void hycopLangevinParse(GROUP *group);
#endif 
