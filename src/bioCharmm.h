#ifndef BIOCHARMM_H
#define BIOCHARMM_H

#include "system.h"
#include "bioCharmmParms.h"
#include "bioCharmmPar.h"

typedef struct res_range_struct
{
    int start;
    int end;
} RESRANGE;

int compare(const void * a, const void * b);

typedef struct charmmlj_struct
{
    double rcut, eps, sigma;
    double shift;
} CharmmLJ_PARMS;


//CharmmLJ_PARMS **CharmmLJ_parms(CHARMMPOT_PARMS *parms);
CharmmLJ_PARMS CharmmLJ_parms(POTENTIAL *potential, LJCH_PARMS *speciesA, LJCH_PARMS *speciesB, double cutoff);

//void CHLennardJones_setShift(CharmmLJ_PARMS *ljparms);
double CHLennardJones_setShift(double sigma, double eps, double rcut);
double CHLennardJones(CharmmLJ_PARMS *ljparms, double r, double *dvdr);


void charmmConvalent(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e);

void printBioEnergies(BIOENERGIES en);

#endif /* BIOCHARMM_H */

