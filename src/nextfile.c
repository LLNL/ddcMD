#include "nextfile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include <wordexp.h>
#include <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "ddc.h"
#include "simulate.h"
#include "system.h"
#include "ddcMalloc.h"
#include "ddcenergy.h"
#include "mpiUtils.h"
#include "state.h"
#include "expandbuffer.h"
#include "gid.h"
#include "collection.h"
#include "pio.h"

PFILE* Pfile_init(const char *filename, const char *mode, MPI_Comm comm);
void   Popen_forRead(PFILE*);
static  pio_long64 next = 0; 

NEXTFILE_PARMS *nextfile_parms(INTEGRATOR*integrator)
{
   NEXTFILE_PARMS *parms = ddcMalloc(sizeof(NEXTFILE_PARMS));
   SIMULATE *simulate = (SIMULATE *)integrator->parent;
   SYSTEM *sys = simulate->system; 
   COLLECTION *collection = sys->collection; 
   
   char **files=NULL;
   int nfiles = object_getv((OBJECT *)integrator, "files", (void *)&files, STRING,IGNORE_IF_NOT_FOUND);
   parms->nfiles = nfiles;
   parms->files = files;
   parms->loopindex = 0;

   collection->size = 0; 
   parms->pfile = Popen(parms->files[0], "r", COMM_LOCAL);
   next = parms->pfile->next; 
   collection->dynamicInfoFromFileHeader=1; 
   collection_read(collection,parms->pfile); 
   Pclose(parms->pfile);
   return parms; 

}
void nextfile(DDC*ddc, SIMULATE*simulate, NEXTFILE_PARMS*p)
{
   char **list = NULL; 
   SYSTEM *sys = simulate->system;
   COLLECTION *collection = sys->collection;
   BOX_STRUCT *box = sys->box;

   collection->size=0; 
   p->pfile = Pfile_init(p->files[0], "r", COMM_LOCAL);
   p->pfile->start = next; 
   Popen_forRead(p->pfile);
   next = p->pfile->next;  
   collection_read(collection,p->pfile); 
   Pclose(p->pfile);

   ddcenergy(ddc, sys, 1);
   kinetic_terms(sys, 1);
   eval_energyInfo(sys);
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
