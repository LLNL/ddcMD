#ifndef DATASUBSET_H
#define DATASUBSET_H
#include <stdio.h>
#include "collection.h"
#include "analysis.h"

typedef struct field_st
{
   int type; 
   char*name;
   char*units;
   char *format; 
   double convert; 
   double accumulatedValue;
   double  (*function)(COLLECTION *, int i);
   double nParticleSamples; 

} FIELD;
typedef struct dataSubset_st 
{
   int nFields; 
   FIELD *fields; 
   char *filename; 
   FILE *file; 
   int nEvals; 
} DATASUBSET_PARMS;

DATASUBSET_PARMS *dataSubset_parms(ANALYSIS *analysis);
void dataSubset_eval(ANALYSIS* analysis);
void dataSubset_output(ANALYSIS* analysis);
void dataSubset_startup(ANALYSIS* analysis);
#endif
