#ifndef PAIR_H_
#define PAIR_H_

#include "potential.h"
#include "neighbor.h"
#include "system.h"
#include "energyInfo.h"

typedef struct pair_parms_str
{
        POTENTIAL *potential; 
        unsigned nspecies; 
        char **species_list; 
        RCUT_TYPE *rcut;
        double rmax;
        double(*fcn)(void *,double, double *); 
        double *rCutoff;
        void **parms; 
	unsigned *type;
	double *q;
} PAIR_PARMS;

typedef struct lj_struct 
{ 
        double rcut,eps,sigma;
   double shift;
}  LJ_PARMS;
void *table_parms(POTENTIAL *potential);

void pair(SYSTEM *sys, PAIR_PARMS *parms, ETYPE *e);
LJ_PARMS **LennardJones_parms(PAIR_PARMS *parms);

#endif //PAIR_H
