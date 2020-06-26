#include <stdio.h>
#include <assert.h>
#include <mpi.h>

#include "object.h"
#include "simulate.h"
#include "three_algebra.h"
#include "ddcMalloc.h"
#include "species.h"
#include "group.h"
#include "units.h"
#include "state.h"
#include "format.h"

int getRank(int);

typedef struct vcmWrite_parms_st
{
   FILE* file;
} VCMWRITE_PARMS;


VCMWRITE_PARMS* vcmWrite_parms(ANALYSIS* analysis)
{
   VCMWRITE_PARMS* parms = ddcMalloc(sizeof(VCMWRITE_PARMS)); 
	
   OBJECT* obj = (OBJECT*) analysis;

   char* filename;

   object_get(obj, "filename", &filename, STRING, 1, "vcm.data");

   if (getRank(0) == 0)
   {
      parms->file = fopen(filename, "a");
      assert(parms->file != NULL);
      
      char fmt[16]; sprintf(fmt,"-%%%ds %%14s",loopFormatSize()); 
      fprintf(parms->file, fmt, "#loop", "time(fs)");
      fprintf(parms->file, "%-51s", "     System vx vy vz (Ang/fs)");
      unsigned nGroups= system_getNgroup(NULL);
      for (unsigned ii=0; ii<nGroups; ++ii)
      {
	 GROUP* g = group_by_index(NULL, ii);
	 assert(g != NULL);
	 char temp[51];
	 sprintf(temp, "     Group %s vx vy vz (Ang/fs)", g->name);
	 fprintf(parms->file, "%-51s", temp);
      }
      unsigned nSpecies = system_getNspecies(NULL);
      for (unsigned ii=0; ii<nSpecies; ++ii)
      {
	 SPECIES* s = species_by_index(NULL, ii);
	 char temp[51];
	 sprintf(temp, "     Species %s vx vy vz (Ang/fs)", s->name);
	 fprintf(parms->file, "%-51s", temp);
      }
      fprintf(parms->file, "\n");
      fflush(parms->file);
   }
   
   ddcFree(filename);
   return parms;
}

/** The vcm array contains data in this order:
 *  vmlocal[0] = whole system
 *  vmlocal[1:1+nGroups] = vmlocal per group
 *  vmlocal[1+nGroups:nGroup+nSpecies] = vmlocal per species
*/
void vcmWrite_output(ANALYSIS* analysis)
{
   VCMWRITE_PARMS* parms = (VCMWRITE_PARMS*) analysis->parms;
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   SYSTEM* sys=simulate->system;
   STATE* state = sys->collection->state;
   FILE* file = parms->file;
   unsigned nlocal = sys->nlocal;
   
//__   double* vx = sys->collection->state->vx;
//__   double* vy = sys->collection->state->vy;
//__   double* vz = sys->collection->state->vz;
   SPECIES** species = sys->collection->state->species;
   GROUP** group = sys->collection->state->group;
   
   unsigned nSpecies=system_getNspecies(NULL);
   unsigned nGroups= system_getNgroup(NULL);
   unsigned nv = nSpecies + nGroups + 1;
   THREE_VECTOR* vmlocal = (THREE_VECTOR*)  ddcCalloc(nv, sizeof(THREE_VECTOR));
   THREE_VECTOR* vmsum = (THREE_VECTOR*)  ddcCalloc(nv, sizeof(THREE_VECTOR));
   double* mlocal = (double*)  ddcCalloc(nv, sizeof(double));
   double* msum = (double*)  ddcCalloc(nv, sizeof(double));
   THREE_VECTOR *__v = state_getVPtr(sys->collection->state); 
   
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      unsigned iSpecies = species[ii]->index;
      unsigned iGroup = group[ii]->index;
      double mass = ((ATOMTYPE_PARMS *)(state->species[ii]->parm))->mass;
      THREE_VECTOR vi;
//__      VSET(vi, vx[ii], vy[ii], vz[ii]);
      vi = __v[ii]; 
      VSCALE(vi, mass);
      VECACUM(vmlocal[0], vi);
      VECACUM(vmlocal[1+iGroup], vi);
      VECACUM(vmlocal[1+nGroups+iSpecies], vi);
      mlocal[0] += mass;
      mlocal[1+iGroup] += mass;
      mlocal[1+nGroups+iSpecies] += mass;
   }
   
   MPI_Reduce(vmlocal, vmsum, 3*nv, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
   MPI_Reduce(mlocal, msum, nv, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);

   if (getRank(0) == 0)
   {
      double time = units_convert(1.0, NULL, "fs") * simulate->time;
      SIGNED64 loop = simulate->loop;
      fprintf(file, loopFormat(), loop);
      fprintf(file, " %16.6f", time);
      double velocity_convert = units_convert(1.0, NULL, "Ang/fs");
      for (unsigned ii=0; ii<nv; ++ii)
      {
	 if (msum[ii] > 0.0)
	    VSCALE(vmsum[ii], 1/msum[ii]);
	 VSCALE(vmsum[ii], velocity_convert);
	 fprintf(file, " %16.6e %16.6e %16.6e", vmsum[ii].x, vmsum[ii].y, vmsum[ii].z);
      }
      fprintf(file, "\n");
      fflush(file);
   }
   ddcFree(vmlocal);
   ddcFree(vmsum);
   ddcFree(mlocal);
   ddcFree(msum);
}

void vcmWrite_eval(ANALYSIS* analysis)
{
}

void vcmWrite_close(ANALYSIS* analysis)
{
   VCMWRITE_PARMS* parms = (VCMWRITE_PARMS*) analysis->parms;
   fclose(parms->file);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
