#ifndef SYSTEM_H
#define SYSTEM_H
#include <stdio.h>
#include "energyInfo.h"
#include "molecule.h" 
#include "species.h" 
#include "group.h" 
#include "potential.h" 
#include "collection.h" 
#include "neighbor.h"
#include "box.h"
#include "gid.h" 
#include "energyInfo.h"
#include "random.h"
enum SYSTEM_CLASS { NORMAL };
enum PRATE_MASK { pPotentialEnergyMask=1, pVirialMask=2, pSionMask=4 };
enum SYSTEM_ENUM  { RETURN_CODE, SYSTEMNAME, SYSTEMTYPE, SYSTEMITYPE, ENERGY, PRESSURE, NLOCAL, NGLOBAL, CURRENT_SYSTEM, SYSTEM_LOOP, SYSTEM_TIME};
typedef struct system_normal_parms_st
{
	void *p;
} SYSTEM_NORMAL_PARMS;

typedef struct system_st
{
	char *name;		/* name of the system */
	char *objclass;
	char *value;
	char *type;		/* type label */
	void  *parent; 
	enum SYSTEM_CLASS itype;	/* integer label for type */
	MOLECULECLASS *moleculeClass;
	SPECIES **species;
	GROUP **group;		
	RANDOM* random;  
	POTENTIAL **potential;
   enum NEIGHBORTABLETYPE neighborTableType; 
	BOX_STRUCT *box;	
	COLLECTION *collection;
   ENERGYINFO *energyInfoObj; 
	NBR *neighbor;
	SIGNED64 loop;
	double time;
	int nspecies;		/* number of the species in the system */
	int ngroup;		   /* number of the groups in the system */
	int npotential;	/* number of potentials in the system */
   gid_type nglobal;
   unsigned nlocal;
   unsigned nremote, nion;	/* number of the objects in the system */
   unsigned nConstraints;
   double p;		/* pressure */
   double t;		/* temperature */
   double v;		/* volume */
   double energy;		/* energy */
   ETYPE energyInfo; 
   int changed;
   int status;
   int return_code;
   int nthreads;
   int nIdleThreads;
   int pCalculate;
   int *pPotentialEnergyRate;
   int *pVirialRate;
   int *pSionRate;
   void (*kineticFunc) (void *sys, int flag); 
   void *parms;		/* system dependent parameters */
} SYSTEM;

#ifdef __cplusplus
extern "C" { 
#endif


   void system_write(char *string);
   SYSTEM *system_init(void *parent,char *name);
   void *system_element(void **ptr, int n, char *name);

   SIGNED64 system_getLoop(SYSTEM*system);
   double system_getTime(SYSTEM*system) ;
   gid_type system_getNglobal(SYSTEM*system) ;
   unsigned system_getNlocal(SYSTEM*system) ;
   int system_getNspecies(SYSTEM *system);
   int system_getNgroup(SYSTEM *system);
   STATE *system_getState(SYSTEM *system);
   SPECIES **system_getSpecies(SYSTEM *system);
   SYSTEM *system_getSystem(SYSTEM *system);
   GROUP **system_getGroup(SYSTEM *system) ;

   int system_get(SYSTEM*system, int get, void *ptr);
   int system_put(SYSTEM*system, int put, void *ptr);

   void  system_putNlocal(SYSTEM*system,unsigned nLocal) ;
   RANDOM* system_getRandom(SYSTEM* system);
   void system_getParticlesEnds(SYSTEM *system, THREE_VECTOR* rminglobal, THREE_VECTOR* rmaxglobal);
   void zeroParticle(SYSTEM *sys); 
   void zeroEType(ETYPE *energyInfo); 
   unsigned system_pCalculate(SYSTEM *sys); 
   enum POT_COMM_MODE system_getCommMode(SYSTEM* sys);


#ifdef __cplusplus
}
#endif



#endif 

/* Local Variables: */
/* tab-width: 3 */
/* End: */
