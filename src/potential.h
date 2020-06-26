#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <stdio.h>
#include "neighbor.h"

enum POTENTIAL_CLASS { NO_POTENTIAL=-1,
                       ZEROPOTENTIAL,
                       MGPT,
                       EAM,
                       EAM1PASS,
                       EAM2PASS,
                       EAM_OPT,
                       EAM_ONEPASS,
                       PAIR,
                       CHARMM,
                       MARTINI,
                       RESTRAINT,
                       EWALD,
                       PLASMA,
                       ORDERSH,
                       ONEBODY,
                       REFLECT,
                       PAIRENERGY,
                       MEAM,
                       MIRRORSYM,
                       LOCALYUKAWA,
                       HYCOP,
                       FMM,  //! Fast Multipole Method.
                       GPU_PAIR
                       };

enum POT_COMM_MODE {POT_ONESIDED, POT_TWOSIDED};

/** set call_fsumX to one if the potential computes "virtual bond
 * forces".  In this case the generic fsumX code can use the virtual
 * bond forces to compute the actual forces on the atoms.  If the
 * potential takes care of forces on the atoms directly then set
 * call_fsumX = 0.*/


typedef struct potential_st
{
   char *name;		/* potential name */
   char *objclass;
   char *value;
   char *type;		/* model */
   void *parent; 
   enum POTENTIAL_CLASS itype;	/* integer label for type */
   void (*eval_potential) (void *sys, void *parms, void *e);	/* pointer to potential function */
   void (*write_dynamics) (void *potential,FILE *file);
   RCUT_TYPE* (*getCutoffs) (void* sys, void* parms, int* nCutoffs);
   enum NEIGHBORTABLETYPE neighborTableType; 
   int call_fsumX;
   int use_gpu_list;
	enum POT_COMM_MODE commMode;
   void *parms;		/* pointer to  parameters for potential function */
} POTENTIAL;


POTENTIAL *potential_init(void *parent,char *name);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
