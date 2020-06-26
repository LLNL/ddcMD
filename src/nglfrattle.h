#ifndef NGLFRATTLE_H
#define NGLFRATTLE_H
#include "ddc.h"
#include "simulate.h"
#include "species.h"
#include "bioCharmm.h"

typedef struct nglfrattle_parms_st
{
    THREE_VECTOR *acc;
    CHARMMPOT_PARMS *charmmpot_parms;
    void (*moveFunc)(SIMULATE*simulate, STATE* state, GID_ORDER* gidOrder, int nTotal,
        char * name, RESRANGE resRange, RESI_CONN* resiConn);
} NGLFRATTLE_PARMS;

NGLFRATTLE_PARMS *nglfrattle_parms(INTEGRATOR*integrator);
void nglfrattle(DDC *ddc, SIMULATE *simulate, NGLFRATTLE_PARMS *p);
NGLFRATTLE_PARMS *nglfrattleTest_parms(INTEGRATOR*integrator);
void nglfrattleTest(DDC *ddc, SIMULATE *simulate, NGLFRATTLE_PARMS *p);
NGLFRATTLE_PARMS *nglfrattleError_parms(INTEGRATOR*integrator);
void nglfrattleError(DDC *ddc, SIMULATE *simulate, NGLFRATTLE_PARMS *p);

typedef struct ref_str
{
    THREE_VECTOR r, v;
} REF;
REF *vaf_v0(void);



void scalePositionsByBoxChange_sk(BOX_STRUCT *box, double time, double *rx, double *ry, double *rz, unsigned nlocal);
void updateStateAliases_sk(SYSTEM *sys, unsigned *nlocal, unsigned *nion, double **rx, double **ry, double **rz, double **vx, double **vy, double **vz, double **fx, double **fy, double **fz, SPECIES ***species, gid_type **label);
int compareSoloGid(const void *v1, const void *v2);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
