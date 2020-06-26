#include "simulate.h"
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <wordexp.h>
#include <inttypes.h>
#include <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "pio.h"
#include "io.h"
#include "system.h"
#include "utilities.h"
#include "printinfo.h"
#include "eq.h"
#include "box.h"
#include "ddcMalloc.h"
#include "transform.h"
#include "units.h"
#include "codata.h"
#include "heap.h"
#include "auxNeighbor.h"
#include "routineManager.h"
#include "ompUtils.h"
#include "mpiUtils.h"
#include "gid.h"
#include "format.h"

static SIMULATE *simulate_last;

int getRank(int);

ANALYSIS  *analysis_init(SIMULATE *simulate,char *name);

int atomsdirParse(char *root, char **atomsdir)
{
	char *dirname;
	wordexp_t buf;
	int i, id, ndir, Index;
	id = getRank(0);
	*atomsdir = NULL;
	_wordexp(root, &buf, 0);
	ndir = buf.we_wordc;
	Index = id%ndir;
	dirname = buf.we_wordv[Index];
	i = strlen(dirname);
	if (dirname[i - 1] == '/') dirname[i - 1] = '\0';
	if (*atomsdir != NULL) ddcFree(*atomsdir);
	*atomsdir = strdup(dirname);
	i = strlen(root);
	if (root[i - 1] == '/') root[i - 1] = '\0';
	_wordfree(&buf);
	return Index;
}


/** Initialization Order:
 *
 *  The SIMULATE object owns both the SYSTEM and the DDC and one of them
 *  must be initialized before the other.  The trouble is that there are
 *  complex dependencies between the DDC and the SYSTEM.  In particular,
 *  the coordinates of the domain centers in the DDC cannot be correctly
 *  initialized without the h-matrix in the BOX that is owned by the
 *  SYSTEM.  This would indicate that the SYSTEM should be initialized
 *  first.  However, the SYSTEM owns the POTENTIAL(s) including possibly
 *  a PPPM object that needs to know about the domain centers in the DDC
 *  in order to build comm routings.  Hence it would be useful to
 *  initialize the DDC before the SYSTEM.
 *
 *  To resolve this problem we now (June 2011) specify the following:
 *
 *  ****  The SYSTEM will be initialized *before* the DDC *****
 *
 *  All code under the call tree of system_init *must* be independent of
 *  the DDC since it will not be initialized.
 *
 *  Rationale:
 *
 *  1) The construction of the h-matrix is non-trivial.  We cannot
 *  assume that the h-matrix is given in the object database.  For
 *  example, in the case of a COLLECTION constructed from a lattice
 *  specification no h-matrix is needed in the input.  If one is
 *  specified it is ignored and overwritten by the lattice vectors times
 *  the lattice size.  Hence, there is no trivial way to obtain a
 *  correct h-matrix other than calling system_init.  This means there
 *  is no reliable method to properly initialize the domain centers
 *  until after the SYSTEM is initialized.
 *
 *  2) Any object or structure in the SYSTEM that needs DDC information
 *  can defer the DDC dependent operations (such as the creation of comm
 *  tables) until such time as those operations are necessary.  It must
 *  be possible to do this since the DDC can't be correctly built during
 *  system_init.
 *
 *  3) If there is any operation in the system_init call tree that needs
 *  DDC info, and cannot be deferred, then we have a serious design
 *  problem that will require further rethinking of code structure.  Our
 *  current philosophy is that the physical system needs no information
 *  about the parallelization.
 */
SIMULATE *simulate_init(void *parent, char *name, MPI_Comm comm)
{
   if (name == NULL && object_exists("simulate", "SIMULATE")) name = "simulate"; 
   if (name == NULL && object_exists("test", "SIMULATE")) name = "test"; 
   if (name == NULL)  error_action("No name provided for  SIMULATE object", ERROR_IN("simulate_init", ABORT));
   if (!object_exists(name, "SIMULATE"))  error_action("A SIMULATE object with name ",name,", can not be found", ERROR_IN("simulate_init", ABORT));
   timestamp("Task 0 started simulate_init");
	char *type, *systemname, *printinfo, *string, *dataBrokerName, *acceleratorName;
	char **names = NULL;
	static SIMULATE *simulate;
	int i, deltaloop;


	simulate = (SIMULATE *) object_initialize(name, "SIMULATE", sizeof(SIMULATE));
	simulate->parent = parent; 
   simulate->comm = comm; 
	simulate_last = simulate; 
   simulate->system=NULL; 

	simulate->ddc = NULL;
	simulate->ensemble = NULL;
	simulate->printinfo = NULL;
	simulate->analysis = NULL;
	simulate->transform = NULL;
	simulate->integrator = NULL;
	simulate->atomsdir = NULL;
	simulate->atomsdirOrginal = NULL;
	simulate->snapshotdir = NULL;
	simulate->Veq = NULL;
   simulate->Peq = NULL;
   simulate->Teq = NULL;
   simulate->dataBroker=NULL;
   simulate->accelerator=NULL;
   simulate->dynamicWrite.n=0; 
   simulate->dynamicWrite.info=ddcMalloc(sizeof(DYNAMICRESTARTINFO)); 

   OBJECT* obj = (OBJECT*) simulate;
   object_get(obj, "heap", &string, STRING, 1, "heap");
   heap_init(simulate,string);  

   //Before the call was:
   //  callRoutine(NULL,slave_heap_init,sizeof(string)+1,string); 
   // I think what is actually meant is:
   callRoutine(NULL,slave_heap_init,strlen(string)+1,string); 

   ddcFree(string); 

   object_get(obj, "type", &(type), STRING, 1, "MD");
   object_get(obj, "loop", &simulate->loop, U64, 1, "0");
   object_get(obj, "maxloop", &simulate->maxloop, U64, 1, "0");
   object_get(obj, "deltaloop", &deltaloop, INT, 1, "-1");
   object_get(obj, "printrate", &simulate->printrate, INT, 1, "5");
   object_get(obj, "snapshotrate", &simulate->snapshotrate, INT, 1, "100");
   object_get(obj, "checkpointrate", &simulate->checkpointrate, INT, 1, "1000");
   object_get(obj, "minrate", &simulate->minrate, INT, 1, "1");
   object_get(obj, "snapshotRootDir", &simulate->atomsdirOrginal, STRING, 1, "./");
   object_get(obj, "time", &simulate->time, WITH_UNITS, 1, "0.0","t",NULL);
   object_get(obj, "dt", &simulate->dt, WITH_UNITS, 1, "1.0","t",NULL);
   object_get(obj, "system", &systemname, STRING, 1, NULL);
   object_get(obj, "printinfo", &printinfo, STRING, 1, "printinfo");
   object_get(obj, "run_id", &simulate->run_id, INT, 1, "0");
   object_get(obj, "nfiles", &simulate->nfiles, INT, 1, "0");
   object_get(obj, "checkpointmode", &string, STRING, 1, "ASCII");
   object_get(obj, "nthreads",&simulate->nthreads, INT, 1, "1"); 
   object_get(obj, "databroker", &(dataBrokerName), STRING, 1, "NotSet");
   object_get(obj, "accelerator", &(acceleratorName), STRING, 1, "NoAccelerator");

   //Initialize Accelerator - GPU
   if(strcmp(acceleratorName,"NoAccelerator") != 0)  simulate->accelerator=accelerator_init(simulate, acceleratorName);        

   simulate->atomsdirIndex = atomsdirParse(simulate->atomsdirOrginal, &simulate->atomsdir);
   int nameLength =strlen(simulate->atomsdir)+20;
   simulate->snapshotdir=ddcMalloc(nameLength);

   char *gidFmt=NULL; 
   int nLoopDigits=8; 
   object_get(obj, "gidFormat", &(gidFmt), STRING, 1, "decimal");
   object_get(obj, "nLoopDigits", &(nLoopDigits), INT, 1, "8");
   gidFormatInit(gidFmt); 
   loopFormatInit(nLoopDigits); 

#ifndef WITH_OMP
   simulate->nthreads = 1;
#endif
   setOmpNumThreads(simulate->nthreads);
   if (strcasecmp(string, "BINARY") == 0) simulate->checkpointmode = BINARY_CHECKPOINT;
   else simulate->checkpointmode = ASCII_CHECKPOINT;
   ddcFree(string);
   object_get(obj, "checkpointprecision", &string, STRING, 1, "FULL");
   if (strcasecmp(string, "FULL") == 0)
      simulate->checkpointprecision = FULL_CHECKPOINT;
   else
      simulate->checkpointprecision = BRIEF_CHECKPOINT;
   ddcFree(string);
   if (getRank(0) == 0 )
   {
      FILE *file = fopen("hpm.data","a"); 
      char *time_str = timestamp_string();
      fprintf(file,"\n******* HPM paragraph started at  %s with loop=%"PRIu64"******\n",time_str,simulate->loop); 
      fclose(file); 
   }


   if (getRank(0) == 0)
      printf("SIMULATE: nthreads = %d\n",simulate->nthreads);
   if (simulate->run_id==0) simulate->run_id = (int) time(NULL);
   if (simulate->nfiles > 0)
   {
      Pio_setNumWriteFiles(simulate->nfiles);
      callRoutine(NULL, slave_Pio_setNumWriteFiles, sizeof(int), &(simulate->nfiles));
   }
   simulate->name = ddcCalloc(strlen(name) + 1, sizeof(char));
   simulate->type = ddcCalloc(strlen(type) + 1, sizeof(char));
   strcpy(simulate->name, name);
   if (strcmp(type, "MD") == 0)        simulate->itype = MD;
   if (strcmp(type, "TRANSFORM") == 0) simulate->itype = SIM_TRANSFORM;
   ddcFree(type);
   object_get(obj, "ensemble", &simulate->ensemble, STRING, 1, "NVE");
   if (strcmp(simulate->ensemble, "NPT") == 0)
   {
      char** list;
      object_get(obj, "Peq", &string, LITERAL, 1, NULL);
      simulate->Peq = eq_parse(string,"m/t^2/l","t");
      if (simulate->Peq->function) (simulate->Peq->function) (simulate->time, (void *)simulate->Peq);
      int nelements = object_getv(obj, "Veq", (void *)&list, STRING,IGNORE_IF_NOT_FOUND);
      if (list != NULL) simulate->Veq = eq_parse(list[0],"l^3","t");
      for (int ii=0; ii<nelements; ++ii)
         ddcFree(list[ii]);
      ddcFree(list);
   }
   auxNeighbor_init(simulate,"auxNeighbor");
   simulate->system = system_init(simulate,systemname);
   ddcFree(systemname);
   simulate->printinfo = printinfo_init(simulate,printinfo);
   ddcFree(printinfo);
   object_get(obj, "integrator", &string, STRING, 1, NULL);
   simulate->integrator = integrator_init(simulate,string);
   ddcFree(string);
   if (deltaloop > -1 ) simulate->maxloop = MIN(simulate->maxloop,simulate->loop+deltaloop);

   object_get(obj, "ddc", &string, STRING, 1, "ddc");
   simulate->ddc=ddc_init(simulate, string);
   ddcFree(string); 
   // This should be part of ddc_init, but DDC doesn't know about SYSTEM.
   simulate->ddc->centersAffineUpdateIndex = box_newAffineUpdateIndex(simulate->system->box);
   // look at the comm mode requested by the potential(s) and pass it on to the DDC.
   if (system_getCommMode(simulate->system) == POT_TWOSIDED)
      simulate->ddc->comm_mode = TWOSIDED;
   else
      simulate->ddc->comm_mode = ONESIDED;

   simulate->dataBroker = NULL;
   if(strcmp(dataBrokerName,"NotSet") !=0)
      simulate->dataBroker = dataBroker_init(simulate, dataBrokerName);
   ddcFree(dataBrokerName);
   dbUpdateInitialParms(simulate->dataBroker, simulate->system);

   timestamp("Finished ddc_init");


   simulate->nanalysis = object_getv(obj, "analysis", (void *)&names, STRING,IGNORE_IF_NOT_FOUND);
   if (simulate->nanalysis > 0)
   {
      simulate->analysis = ddcCalloc(simulate->nanalysis, sizeof(ANALYSIS *));
      for (i = 0; i < simulate->nanalysis; i++)
      {
         simulate->analysis[i] = analysis_init(simulate,names[i]);
         ddcFree(names[i]);
      }
      ddcFree(names);
   }
   simulate->ntransform = object_getv(obj, "transform", (void *)&names, STRING,IGNORE_IF_NOT_FOUND);
   if (simulate->ntransform > 0)
   {
      simulate->transform = ddcCalloc(simulate->ntransform, sizeof(TRANSFORM *));
      for (i = 0; i < simulate->ntransform; i++)
      {
         simulate->transform[i] = transform_init(simulate,names[i]);
         ddcFree(names[i]);
      }
      ddcFree(names);
   }

   //simulate->atomsdirIndex = atomsdirParse(simulate->atomsdirOrginal, &simulate->atomsdir);
   //int nameLength =strlen(simulate->atomsdir)+20;
   //simulate->snapshotdir=ddcMalloc(nameLength);
   SYSTEM *sys = simulate->system; 
   STATE* state = sys->collection->state;
   state->nlocal = sys->nlocal;
   state->nion = sys->nion;

   timestampBarrier("All tasks finished simulate_init",COMM_LOCAL);
   return simulate;
}

int simulate_getCheckpointRate(SIMULATE *simulate) { if (simulate == NULL) simulate=simulate_last; return simulate->checkpointrate;}
int64_t simulate_getMaxLoop(SIMULATE *simulate) { if (simulate == NULL) simulate=simulate_last; return simulate->maxloop;}
int64_t simulate_getLoop(SIMULATE *simulate) { if (simulate == NULL) simulate=simulate_last; return simulate->loop;}
double simulate_getTime(SIMULATE *simulate) { if (simulate==NULL) simulate=simulate_last;return simulate->time;}
int simulate_getNthreads(SIMULATE *simulate) { if (simulate==NULL) simulate=simulate_last;return simulate->nthreads;}
SIMULATE *simulate_getSimulate(SIMULATE *simulate) { if (simulate==NULL) simulate=simulate_last;return simulate;}
int simulate_get(SIMULATE*simulate, int get, void **ptr)
{
   if (simulate == NULL) simulate = simulate_last; 
   switch (get)
   {
      case SIMULATENAME:
         *ptr = (void *)simulate->name;
         return 1;
      case SIMULATEITYPE:
         *((enum SIMULATE_CLASS *)ptr) = simulate->itype;
         return 1;
      case SIMULATETYPE:
         *ptr = (void *)simulate->type;
         return 1;
      case SIMULATE_LAST:
         *ptr = (void *)simulate_last; 
         return 1; 
      default:
         return 0;
   }
   return 0;
}

int simulate_put(SIMULATE*simulate, int put, void *ptr)
{
   if (simulate == NULL) simulate = simulate_last; 
   switch (put)
   {
      case SIMULATE_LOOP:
         simulate->loop = *(SIGNED64*)ptr; 
         system_put(simulate->system,SYSTEM_LOOP,(void *)&simulate->loop); 
         return 1 ; 
      case SIMULATE_TIME:
         simulate->time = *(double*)ptr; 
         system_put(simulate->system,SYSTEM_TIME,(void *)&simulate->time); 
         return 1 ; 
      default:
         return 0;
   }
}
void registerDynamicWrite(void *parms, void (*fcn)(void*, FILE *file))
{
   SIMULATE *simulate = simulate_getSimulate(NULL);
   int n = simulate->dynamicWrite.n; 
   simulate->dynamicWrite.info =
      ddcRealloc(simulate->dynamicWrite.info,(n+1)*sizeof(DYNAMICRESTARTINFO)); 
   simulate->dynamicWrite.n++; 
   simulate->dynamicWrite.info[n].parms=parms; 
   simulate->dynamicWrite.info[n].fcn=fcn; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
