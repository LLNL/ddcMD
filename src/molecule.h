#ifndef MOLECULE_H
#define MOLECULE_H
#include "gid.h"
#include "species.h"
#include "state.h"
enum MOLECULE_CLASS { NORMALMOLECULE, GIDORDERED };
typedef struct moleculeType_st
{
	char *name;		/* molecule name */
	char *objclass;
	char *value;
	char *type;		/* model */
	void *parent; 
	enum MOLECULE_CLASS itype;	/* integer label for type */
	int index;
	int nSpecies;
   SPECIES **species;
   SPECIES *ownershipSpecies; 
   int ownershipSpeciesOffset; 
	int nMembers; 
	void *parm;		/* pointer to  parameter */
} MOLECULETYPE;
typedef struct  molecule_st
{
  unsigned gid; 
  MOLECULETYPE *type; 
  int *list; 
} MOLECULE;

typedef struct moleculeClass 
{
	char *name;		/* molecule name */
	char *objclass;
	char *value;
	char *type;		/* model */
	void *parent; 
   int nMoleculeTypes;
   MOLECULETYPE **moleculeTypes; 
   int nMoleculesGlobal; 
   int nMolecules; // number of multispecies molecule. 
   int nMoleculesSpecies; 
   MOLECULE *molecule;
   int nSpecies; 
   SPECIES **species; 
   SPECIES **speciesBuffer;
   int *speciesIndexToMoleculeIndex;
   int *speciesToMoleculeOffset;
   int *list;
   
   int orderType; 
   int *orderList;    // 
} MOLECULECLASS;

MOLECULECLASS *moleculeClassInit(void *parent, char *name);
void moleculeScanState(MOLECULECLASS *moleculeClass, STATE *state);
void moleculeDelta(MOLECULECLASS *moleculeClass, STATE *state, THREE_VECTOR *delta);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
