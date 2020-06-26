#ifndef LANGEVIN_H
#define LANGEVIN_H
#include "gid.h"
#include "eq.h"
#include "three_algebra.h"
#include "lcg64.h"
#include "random.h"
#include "group.h"
#include "potential.h"
enum LANGEVIN_TYPE { NORMALLANGEVIN, SKUPSKYLANGEVIN, LOCALYUKAWALANGEVIN, HYCOPLANGEVIN };
enum LANGEVIN_TEMP { UNIFORM, EXPLICIT_TIME, GLOBAL_ENERGY};
#define LANGEVIN_PARMS_COMMON    int type, stridekBT, strideTau; double *tau; double *kBT; THREE_VECTOR vcm; RANDOM *random; 
typedef struct langevin_parms_st
{
   LANGEVIN_PARMS_COMMON 

	EQTARGET *Teq_vs_time;
	double Teq; 
   int Teq_dynamics;
   int total_energy_set;
	double total_energy;
   double Cp; 
//   THREE_VECTOR *g; 

} LANGEVIN_PARMS;

/*
typedef struct localYukawaLangevin_parms_st
{
   LANGEVIN_PARMS_COMMON; 

	EQTARGET *Teq_vs_time;
	double Teq; 
   char *potentialName; 
   POTENTIAL *potential; 
} LOCALYUKAWALANGEVIN_PARMS;

*/
void langevin_write_dynamics(GROUP*g, FILE *file);
double *langevin_allocateTau(GROUP *gp ,int n);
double langevin_getTemperature(LANGEVIN_PARMS *parms, double time);
double langevin_getTauType(GROUP *g);
double *langevin_getTau(GROUP *g);
#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
