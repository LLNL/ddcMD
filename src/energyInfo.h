#ifndef ETYPE_H
#define ETYPE_H
#include "three_algebra.h"

enum ENERGYINFO_CLASS { ENERGYINFO_NORMAL} ;
typedef struct energyInfo_st
{
	char *name;		/* name of the system */
	char *objclass;
	char *value;
	char *type;		/* type label */
	void  *parent; 
	enum ENERGYINFO_CLASS itype;	/* integer label for type */
   unsigned *groupEvalRate;
   unsigned *globalEvalRate;
} ENERGYINFO;

typedef struct energyInfoType_st
{
   int called;
	double number;
   double mass; 
	double rk_local, eion_local, pion_local;
   double rk, eion, pion;
   double eBath; 
   double e0_at_V; 
   double pdiaV;
   double vol_per_atom; 

	THREE_VECTOR vcm;
   THREE_VECTOR f;
   THREE_VECTOR thermal_flux; 
	THREE_SMATRIX virial; 
	THREE_SMATRIX sion; 
	THREE_SMATRIX sion_local; 
   THREE_SMATRIX tion;
   THREE_VECTOR  *virialCorrection; 

   double temperature;
} ETYPE;
struct system_st; 
void eval_energyInfo(struct system_st *sys);
ENERGYINFO *energyInfo_init(void *parent, char *name);
#endif 



/* Local Variables: */
/* tab-width: 3 */
/* End: */
