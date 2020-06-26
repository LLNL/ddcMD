
/*
   enum PRINTINFO_ENUM { PLENGTH,PMASS,PTIME,PCURRENT,PTEMPERATURE,PAMOUNT,PLUMINOUS_INTENSITY,PENERGY,PPRESSURE,PVOLUME,PFORCE,PENERGYFLUX,NUNITS}; 
   static char *unitnames[] =         { "LENGTH","MASS","TIME","CURRENT","TEMPERATURE","AMOUNT","LUMINOUS_INTENSITY","ENERGY","PRESSURE","VOLUME","FORCE","ENERGYFLUX"}; 
   static char *unitDefaultValue[] = { "Ang",   "amu", "fs",  "e/fs",   "K",          "1",     "1",                  "eV",    "GPa",     "Bohr^3","eV/Ang","ueV/Ang^2/fs"}; 
   for (int i=0;i<NUNITS;i++) 
   {
   printinfo->units[i].name=unitnames[i];
   object_get((OBJECT*) printinfo,unitnames[i],&printinfo->units[i].value,STRING,1,unitDefaultValue[i]);
   }
 */

// accumulate a simple time average of the force on certain particles
// append the values to an output file
//
#include "dataSubset.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

//#include "box.h"
#include "object.h"
#include "pio.h"
#include "io.h"
#include "crc32.h"
#include "ddcMalloc.h"
#include "units.h"
#include "mpiAlgorithm.h"
#include "pinfo.h"
#include "ioUtils.h"
#include "system.h"
#include "three_algebra.h"
#include "mpiUtils.h"
#include "format.h"


int getRank(int);
void dataSubset_clear(ANALYSIS* analysis);

typedef double (* dataSubsetFunctionType)(COLLECTION*,int ) ;


static char *scat512(char *a,char*b)
{ 
   static char s[512];
   snprintf(s,511,"%s(%s)",a,b);
   return s;
}


double particleEnergy(COLLECTION *collection,  int i)
{
   STATE *state= collection->state; 
   double mass = ((ATOMTYPE_PARMS *)(state->species[i]->parm))->mass;
   double vx = state->vx[i];
   double vy = state->vy[i];
   double vz = state->vz[i];
   double u =  state->potentialEnergy[i];
   double e =  0.5*mass*(vx*vx+vy*vy+vz*vz)+u; 
   return e; 
}
double particlePotential(COLLECTION *collection,  int i)
{
   STATE *state= collection->state; 
   double u =  state->potentialEnergy[i];
   return u; 
}
double particleKinetic(COLLECTION *collection,  int i)
{
   STATE * state= collection->state; 
   double mass = ((ATOMTYPE_PARMS *)(state->species[i]->parm))->mass;
   double vx = state->vx[i];
   double vy = state->vy[i];
   double vz = state->vz[i];
   double k =  0.5*mass*(vx*vx+vy*vy+vz*vz); 
   return k;
}

enum fieldList { TIME, NSAMPLES, NPARTICLES, Etotal, Ekinetic, Epotential, Rx, Ry, Rz, Vx, Vy, Vz, Fx, Fy, Fz}; 
static char *_fieldNames[]= { "time", "nSamples","nParticles","Etotal", "Ekinetic", "Epotential", "Rx", "Ry", "Rz", "Vx", "Vy", "Vz", "Fx", "Fy", "Fz",NULL}; 
static char *_fieldUnits[]= { "fs"  ,"1"         , "1"        ,"eV", "eV", "eV", "Ang", "Ang", "Ang", "Ang/fs", "Ang/fs", "Ang/fs", "eV/Ang", "eV/Ang", "eV/Ang",NULL}; 
static dataSubsetFunctionType _fieldFunctions[] = {NULL,NULL,NULL,particleEnergy,particleKinetic,particlePotential}; 
static char *_fieldFormat[] = {" %16.6f"," %16.0f"," %16.0f"," %16.12f"," %16.12f"," %16.12f"," %16.8f"," %16.12f", " %15.12f"};

static int fieldtype(char *name)
{
   int i=0; 
   while (_fieldNames[i] != NULL) 
   {
      if (strcmp(name,_fieldNames[i])==0) return i; 
      i++; 
   }
   return -1; 
}

DATASUBSET_PARMS *dataSubset_parms(ANALYSIS *analysis)
{
   DATASUBSET_PARMS *parms = ddcMalloc(sizeof(DATASUBSET_PARMS)); 


   OBJECT* obj = (OBJECT*) analysis;
   char *names[256]; 
   char *species[256]; 

   int nNames=object_get(obj, "fields", &names ,STRING, 256, "time nSamples nParticles Etotal Ekinetic Epotential");
   int nSpecies=object_get(obj, "species", &species ,STRING, 256, "");
   char *defaultFilename=ddcMalloc(strlen(analysis->name)+5); 
   sprintf(defaultFilename,"%s.data",analysis->name); 
   object_get(obj, "filename", &parms->filename ,STRING, 1,defaultFilename);
   assert(nNames  < 256);
   assert(nSpecies< 256);
   parms->fields = ddcMalloc(nNames*sizeof(*parms->fields)); 
   parms->nFields = nNames; 
   char formatBuffer[64]; 
   for (int i = 0;i<nNames;i++)
   {
      char *name=strtok(names[i]," ("); 
      int type = fieldtype(name); 
      assert(type != -1); 
      char *units=strtok(NULL," )"); 
      if (units == NULL) units = _fieldUnits[i]; 
      char *format=strtok(NULL," %"); 
      if (format == NULL)  format = _fieldFormat[i]; 
      else {sprintf(formatBuffer," %%%s",format); format = formatBuffer;}
      
      parms->fields[i].type = type; 
      parms->fields[i].name = strdup(name); 
      parms->fields[i].units = strdup(units); 
      parms->fields[i].convert = units_convert(1.0,NULL,units); 
      parms->fields[i].format = strdup(format); 
      parms->fields[i].function = _fieldFunctions[i]; 
      parms->fields[i].accumulatedValue = 0.0; 
      parms->fields[i].nParticleSamples=0; 
   }
   parms->nEvals=0;
   return parms;
}

////////////////////////////////////////////////////////////////////////
void dataSubset_eval(ANALYSIS* analysis)
{

   DATASUBSET_PARMS* parms = (DATASUBSET_PARMS *)analysis->parms; 

   SIMULATE* simulate =(SIMULATE *)analysis ->parent; 
   COLLECTION *collection = simulate->system->collection; 
   int nlocal = simulate->system->nlocal;
   parms->nEvals++; 
   for (int j=0;j<parms->nFields;j++)
   {
      FIELD *field = parms->fields+j; 

      if (field->type == TIME) { field->accumulatedValue    +=  simulate->time; field->nParticleSamples+=1;continue; }
      if (field->type == NSAMPLES) { field->accumulatedValue += 1; field->nParticleSamples+=getSize(0);continue; }
      if (field->type == NPARTICLES) {field->accumulatedValue +=  nlocal; field->nParticleSamples=+1;continue; }
      field->nParticleSamples +=  nlocal;
      for (int i=0;i<nlocal;i++)
      {
            field->accumulatedValue +=  field->function(collection,i); 
      } 
   }
   return ;
}

////////////////////////////////////////////////////////////////////////
void dataSubset_output(ANALYSIS* analysis)
{
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 

   DATASUBSET_PARMS* parms = (DATASUBSET_PARMS *)analysis->parms; 
   
   FILE *file = parms->file; 
   fprintf(file,loopFormat(),simulate->loop); 
   
   for (int j=0;j<parms->nFields;j++)
   {
      FIELD *field = parms->fields+j; 
      double scale = field->convert/field->nParticleSamples; 
      fprintf(file,field->format,scale*field->accumulatedValue); 
   }
   fprintf(file,"\n");  fflush(file); 

   dataSubset_clear(analysis);

}

////////////////////////////////////////////////////////////////////////
void dataSubset_clear(ANALYSIS* analysis)
{
   DATASUBSET_PARMS* parms = (DATASUBSET_PARMS *)analysis->parms;
   for (int i = 0;i<parms->nFields;i++)
   {
      parms->fields[i].accumulatedValue = 0.0; 
      parms->fields[i].nParticleSamples=0; 
   }
   parms->nEvals=0; 
}
void dataSubset_startup(ANALYSIS* analysis)
{
   DATASUBSET_PARMS* parms = (DATASUBSET_PARMS *)analysis->parms;
   parms->file = fopen(parms->filename,"a"); 
   char fmt[8]; sprintf(fmt,"-%%%ds",loopFormatSize()); 
   fprintf(parms->file,fmt,"#loop"); 
   for (int i=0;i<parms->nFields;i++) 
   {
      fprintf(parms->file," %16s",scat512(parms->fields[i].name, parms->fields[i].units));
   }
   fprintf(parms->file,"\n"); 
   dataSubset_clear(analysis);
}

