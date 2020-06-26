#ifndef NUCLEARRX_H
#define NUCLEARRX_H

#include "system.h"
#include "species.h"
#include "nearPair.h"
enum  NUCLEARRX_CLASS { NORMAL_NRX}; 

typedef struct nuclearproducts_st {  int nProducts ; SPECIES **species ; double branchProb,  energy; } nuclearProducts; 
typedef struct CrossSectionParms_st { char *name; SPECIES *rx1 ; SPECIES *rx2; double kqq,reducedMass; int nBranches; nuclearProducts *branches; double bG,A[5],B[5];} CrossSectionParms; 

typedef struct nuclearRx_st
{
   char *name;		/* potential name */
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   enum NUCLEARRX_CLASS itype;	/* integer label for type */
	int nRx; 
	char **names; 
	CrossSectionParms **crossSectionMap;
	FILE *file; 
} NUCLEARRX;


NUCLEARRX *nuclearRx_init(void *parent,char *name);
int nuclearRx_eval(NUCLEARRX *nuclearRx, SYSTEM *sys,NEARPAIR_PARMS *nearPairParms,THREE_VECTOR *R,THREE_VECTOR *V,THREE_VECTOR *F,THREE_VECTOR *A);
void nuclearCrossSection_init(SPECIES *species_n, SPECIES *species_p, SPECIES *species_D, SPECIES *species_T, SPECIES *species_He3, SPECIES *species_He4);
CrossSectionParms *getNuclearCrossSection(SPECIES *sp1, SPECIES *sp2) ;
double lamdaEval(CrossSectionParms *p, double e);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
