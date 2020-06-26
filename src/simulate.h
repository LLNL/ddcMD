#ifndef SIMULATE_H
#define SIMULATE_H
#define  __STDC_FORMAT_MACROS
#include "system.h" 
#include "integrator.h" 
#include "analysis.h"
#include "transform.h"
#include "accelerator.h"
//#include "printinfo.h"
#include "ddc.h"
#include "cdbInterface.h"

enum SIMULATE_CLASS { MD, SIM_TRANSFORM };
enum CHECKPOINT_WRITEMODE {ASCII_CHECKPOINT, BINARY_CHECKPOINT};
enum CHECKPOINT_PRECISION {FULL_CHECKPOINT, BRIEF_CHECKPOINT};
enum SIMULATE_ENUM  { SIMULATENAME, SIMULATETYPE, SIMULATEITYPE, SIMULATE_LAST,SIMULATE_LOOP,SIMULATE_TIME};

typedef  struct dynamicRestartInfo_st
{
   void *parms;
   void (*fcn)(void *,FILE *file); 
} DYNAMICRESTARTINFO ; 
typedef struct dynamicWrite_str
{
   int n; 
   DYNAMICRESTARTINFO *info;
} DYNAMICWRITE; 
typedef struct simulate_st
{
   char *name;		/* name of the system */
   char *objclass;
   char *value;
   char *type;		/* type label */
   void  *parent; 
   MPI_Comm comm;
   DDC *ddc; 
   char *ensemble;
   enum SIMULATE_CLASS itype;	/* integer label for type */
   struct printinfo_st *printinfo; 
   SYSTEM *system;		/* species */
   int nanalysis; 
   ANALYSIS **analysis;
   int ntransform;
   TRANSFORM** transform;
   INTEGRATOR *integrator;
   DATABROKER *dataBroker;
   ACCELERATOR *accelerator;
   int nthreads;
   // start of io parameters
   int run_id; 
   int64_t startLoop;
   int64_t loop; 
   int64_t maxloop; 
   int printrate, snapshotrate, checkpointrate,minrate, atomsdirIndex, nfiles;
   enum CHECKPOINT_WRITEMODE checkpointmode;
   enum CHECKPOINT_PRECISION checkpointprecision;
   char *checkpointUnits;
   char *atomsdir, *atomsdirOrginal, *snapshotdir; 
   DYNAMICWRITE dynamicWrite; 
   // end of io parameters
   double time; 
   double dt;
   EQTARGET *Veq, *Peq, *Teq;
} SIMULATE;

SIMULATE *simulate_init(void *parent, char *name, MPI_Comm comm);
int simulate_getCheckpointRate(SIMULATE *simulate) ;
int64_t simulate_getMaxLoop(SIMULATE *simulate) ;
int64_t simulate_getLoop(SIMULATE *simulate) ;
double simulate_getTime(SIMULATE *simulate) ;
int simulate_getNthreads(SIMULATE *simulate) ;
SIMULATE *simulate_getSimulate(SIMULATE *simulate) ;
int simulate_get(SIMULATE*simulate, int get, void **ptr);
void registerDynamicWrite(void *parms, void (*fcn)(void*, FILE *file));
int simulate_put(SIMULATE*simulate, int put, void *ptr);

#endif 


/* Local Variables: */
/* tab-width: 3 */
/* End: */
