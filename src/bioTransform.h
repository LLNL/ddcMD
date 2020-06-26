#ifndef BIOTRANSFORM_H
#define BIOTRANSFORM_H

#include "bioCharmmTop.h"
#include "cdbInterface.h"

typedef struct transform_str
{
    //Allow molecules to transform?
    //int weightedCharmm;
    //4 possible weights, this holds 3 the other is 0.0
    double * weights;

    //Use files or cdb to update weights?
    int useFiles;
} BIOWEIGHTS;

int charmmTransformInit(BIOWEIGHTS*);
int updateWeights(BIOWEIGHTS*);

double get2Weights(BIOWEIGHTS*, TATOM_PARMS*, TATOM_PARMS*);
double get3Weights(BIOWEIGHTS*, TATOM_PARMS*, TATOM_PARMS*, TATOM_PARMS*);
double get4Weights(BIOWEIGHTS*, TATOM_PARMS*, TATOM_PARMS*, TATOM_PARMS*, TATOM_PARMS*);

#endif 
