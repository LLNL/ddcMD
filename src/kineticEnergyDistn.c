#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>

#include "object.h"
#include "simulate.h"
#include "three_algebra.h"
#include "ddcMalloc.h"
#include "species.h"
#include "group.h"
#include "units.h"
#include "system.h"
#include "species.h"
#include "state.h"
#include "io.h"
#include "format.h"

int getRank(int);

typedef struct bin_st
{
	char *species; 
	int nBins; 
	double minValue,maxValue,delta; 
	double *cnt; 
	double ave; 
	double subCnt,supCnt; 
	double min,max; 
        double cntTotal; 
} BIN; 
typedef struct kineticEnergyDistn_parms_st
{
   int nDistGroups; 
   SPECIES **species; 
   char **distGroups; 
   int *mapS2D; 
   BIN *bin; 
   FILE *file; 
   
} KINETICENERGYDISTN_PARMS;


void kineticEnergyDistn_clear(ANALYSIS* analysis);
KINETICENERGYDISTN_PARMS* kineticEnergyDistn_parms(ANALYSIS* analysis)
{
   KINETICENERGYDISTN_PARMS* parms = ddcMalloc(sizeof(KINETICENERGYDISTN_PARMS)); 
   OBJECT* obj = (OBJECT*) analysis;
   SIMULATE* simulate = (SIMULATE*) analysis->parent; 
   SYSTEM* sys = simulate->system; 


   parms->distGroups=NULL; 
   parms->mapS2D = ddcMalloc(sys->nspecies*sizeof(int)); 
   for (int i=0;i<sys->nspecies;i++) parms->mapS2D[i]=-1; 
   parms->nDistGroups  = object_getv(obj, "distGroups",     (void*)&parms->distGroups, STRING,     IGNORE_IF_NOT_FOUND);
   BIN *bin=parms->bin=ddcMalloc(parms->nDistGroups*sizeof(BIN)); 
   for (int j=0;j<parms->nDistGroups;j++)  
   { 
	OBJECT *objG = object_find(parms->distGroups[j],"BIN"); 
	object_get(objG, "species", &bin[j].species, STRING, 1, "");
	object_get(objG, "emin", &bin[j].minValue, WITH_UNITS, 1, "0", "energy", NULL);
	object_get(objG, "emax", &bin[j].maxValue, WITH_UNITS, 1, "0", "energy", NULL);
	object_get(objG, "nBins", &bin[j].nBins, INT, 1, "1");
	bin[j].cnt = ddcMalloc(bin[j].nBins*sizeof(double)); 
	for (int i = 0;i<bin[j].nBins;i++) bin[j].cnt[i]=0.0; 
	bin[j].delta = (bin[j].maxValue-bin[j].minValue)/bin[j].nBins; 
	bin[j].cntTotal=0.0; 
	bin[j].subCnt=0.0; 
	bin[j].supCnt=0.0; 
	bin[j].ave=0.0; 
	bin[j].min=1e300; 
	bin[j].max=0.0; 
   	for (int i=0;i<sys->nspecies;i++) 
	{
		if (strcmp(parms->bin[j].species,sys->species[i]->name) ==0) 
		{
			assert(parms->mapS2D[i] == -1) ; 
			parms->mapS2D[i] = j; 
			break ; 
		}
	}
   } 
   if (getRank(0) ==0) 
   {
     parms->file=fopen("kinetic.data","a"); 
     fprintf(parms->file,"# loop  time(fs)   \n"); 
   }
   
	

   return parms;
}
static void writeHeader(FILE *file)
{
      fprintf(file, "%-14s %14s %14s\n", "# Energy (eV)", "pdf (1/eV)","cnt");
}
void kineticEnergyDistn_output(ANALYSIS* analysis)
{
   SIMULATE* simulate = (SIMULATE*) analysis->parent; 
   KINETICENERGYDISTN_PARMS* parms = (KINETICENERGYDISTN_PARMS*) analysis->parms;
   CreateSnapshotdir(simulate, NULL);
   for (int i=0;i<parms->nDistGroups;i++) 
   { 
      double buffer[parms->bin[i].nBins+3]; 
      double bufferSum[parms->bin[i].nBins+3]; 
      int j; 
      for (j=0;j<parms->bin[i].nBins;j++) buffer[j] = parms->bin[i].cnt[j]; 
      buffer[j++] = parms->bin[i].cntTotal; 
      buffer[j++] = parms->bin[i].subCnt; 
      buffer[j++] = parms->bin[i].supCnt; 
      buffer[j++] = parms->bin[i].ave; 
      MPI_Reduce(buffer, bufferSum, j, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
      double minBuffer,maxBuffer;
      MPI_Reduce(&parms->bin[i].min, &minBuffer, 1, MPI_DOUBLE, MPI_MIN, 0, COMM_LOCAL);
      MPI_Reduce(&parms->bin[i].max, &maxBuffer, 1, MPI_DOUBLE, MPI_MAX, 0, COMM_LOCAL);
      if (getRank(0) == 0)
      {
         for (j=0;j<parms->bin[i].nBins;j++) parms->bin[i].cnt[j]=bufferSum[j]; 
         parms->bin[i].cntTotal=bufferSum[j++]; 
         parms->bin[i].subCnt=bufferSum[j++]; 
         parms->bin[i].supCnt=bufferSum[j++]; 
         if (parms->bin[i].cntTotal) parms->bin[i].ave=bufferSum[j++]/parms->bin[i].cntTotal; 
         parms->bin[i].min = minBuffer; 
         parms->bin[i].max = maxBuffer; 
      }
   } 
   if (getRank(0) == 0)
   {
      double eC = units_convert(1.0, NULL, "eV");
      for (int i=0;i<parms->nDistGroups;i++) 
      {
         char filename[1024];
         snprintf(filename, 1024,"%s/%s_kDist.data", simulate->snapshotdir,parms->distGroups[i]);
         FILE* file = fopen(filename, "w");
         writeHeader(file); 
         BIN bin = parms->bin[i]; 
         for (int j=0; j<bin.nBins; j++)
         {
            double cnt = bin.cnt[j];
            double cntTotal = bin.cntTotal; 	
            double energy = ((j+0.5)*bin.delta + bin.minValue)*eC; 
            double pdf = cnt/(cntTotal*bin.delta)/eC; 
            fprintf(file,"%e %e %e\n",energy,pdf,cnt); 
         }
         fclose(file);
         fprintf(parms->file,loopFormat(),simulate->loop);
         fprintf(parms->file," %16.6f ", simulate->time);
         fprintf(parms->file,"%12.6f %12.8f %12.8f %4.0f %4.0f %8.0f ", bin.ave*eC,bin.min*eC,bin.max*eC,bin.subCnt,bin.supCnt,bin.cntTotal);
      }
      fprintf(parms->file,"\n"); 
      fflush(parms->file); 
   }
   kineticEnergyDistn_clear(analysis);
}

void kineticEnergyDistn_eval(ANALYSIS* analysis)
{
   SIMULATE* simulate = (SIMULATE*) analysis->parent; 
   KINETICENERGYDISTN_PARMS* parms = (KINETICENERGYDISTN_PARMS*) analysis->parms;
   SYSTEM* sys = simulate->system; 

   double *vx = sys->collection->state->vx;
   double *vy = sys->collection->state->vy;
   double *vz = sys->collection->state->vz;
   SPECIES **species =  sys->collection->state->species;
   unsigned nlocal = sys->nlocal;
   for (unsigned k = 0; k < nlocal; k++)
   {
      int index = species[k]->index; 
      int iDist = parms->mapS2D[index];
      if (iDist == -1) continue; 
      BIN *b = parms->bin+iDist; 
      double mass = ((ATOMTYPE_PARMS *) (species[k]->parm))->mass;
      double v2 = vx[k]*vx[k]+vy[k]*vy[k]+vz[k]*vz[k]; 
      double K = 0.5*mass*(v2);
      b->ave += K; 
      b->cntTotal++; 
      if (K > b->max) b->max = K; 
      if (K < b->min)  b->min = K; 
      if (K <  b->minValue) { b->subCnt+=1.0; continue; }
      if (K >= b->maxValue) { b->supCnt+=1.0; continue; }
      int ibin =  (K-b->minValue)/b->delta; 
      //printf("dist: %d %d %d %d %e %e %e %e\n",k,iDist,ibin,b->nBins,K,b->minValue,b->maxValue,b->delta); 
      assert(ibin >=0 && ibin <b->nBins); 
      b->cnt[ibin]+=1.0; 
   }
}
void kineticEnergyDistn_clear(ANALYSIS* analysis)
{
   KINETICENERGYDISTN_PARMS* parms = (KINETICENERGYDISTN_PARMS*) analysis->parms;
   for (int i=0;i<parms->nDistGroups;i++) 
   {
      BIN *b = parms->bin+i; 
      for (int j=0; j<b->nBins; j++) b->cnt[j]=0.0; 
      b->cntTotal=0.0; 	
      b->subCnt=0.0; 	
      b->supCnt=0.0; 	
      b->min = 1e300; 
      b->max = 0.0; 
      b->ave=0.0; 	
   }
}

void kineticEnergyDistn_close(ANALYSIS* analysis)
{
   //KINETICENERGYDISTN_PARMS* parms = (KINETICENERGYDISTN_PARMS*) analysis->parms;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
