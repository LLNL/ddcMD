#ifndef SPECIES_H
#define SPECIES_H
#include "energyInfo.h"
#include "gid.h"
enum SPECIES_CLASS { ATOMTYPE };
enum SPECIES_ENUMS { NATOMS, NSPECIES, SPECIESTYPE, SPECIESITYPE, SPECIESLIST, NAMELIST, TYPELIST , SPECIESMAXNAMELENGTH, SPECIESMAXTYPELENGTH, SPECIESMAXWRITELENGTH, SPECIESMAXBINARYWRITELENGTH};
typedef struct atomtype_parms_st
{
	double mass, charge;
} ATOMTYPE_PARMS;
typedef struct species_st
{
	char *name;		/* species name */
	char *objclass;
	char *value;
	char *type;		/* model */
	void *parent; 
	enum SPECIES_CLASS itype;	/* integer label for type */
	int index;
	int natoms;
	gid_type nMember; 
	ETYPE  energyInfo; 
	void *parm;		/* pointer to  parameter */
} SPECIES;
SPECIES *species_init(void *parent, char *name);
void species_get(SPECIES*species, int get, void *ptr);
int species_put(SPECIES*species, int put, void **ptr);
SPECIES *species_find(SPECIES ** species, char *name);
SPECIES *species_by_index(SPECIES ** species, int index);
void resetStaticSpeciesVars();
#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
