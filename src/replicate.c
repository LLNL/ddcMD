#include <string.h>
#include "replicate.h"

#include "transform.h"
#include "object.h"
#include "simulate.h"
#include "ddcMalloc.h"
#include "utilities.h"
#include "pio.h"
#include "collection.h"
#include "collection_write.h"
#include "io.h"
#include "box.h"
#include "units.h"
#include "random.h"
#include "mpiUtils.h"
#include "format.h"

static void arrayShift(double* x, unsigned n, double dx);
static void shuffle(gid_type* a, unsigned n);


typedef struct replicate_parms_st
{
   int nx;
   int ny;
   int nz;
   gid_type stride; 
   double gamma; // coefficient of velocity tweak.
   LONG64 seed;
} REPLICATE_PARMS;


void* replicate_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   REPLICATE_PARMS* parms = ddcMalloc(sizeof(REPLICATE_PARMS));
   transform->parms = parms;
   transform->writeNeeded = 0;
   
   int randomize = 0;
   object_get(obj, "nx", &parms->nx, INT, 1, "1");
   object_get(obj, "ny", &parms->ny, INT, 1, "1");
   object_get(obj, "nz", &parms->nz, INT, 1, "1");
   object_get(obj, "stride", &parms->stride, U64, 1, "1");
   object_get(obj, "gamma", &parms->gamma, DOUBLE, 1, "0.0");
   object_get(obj, "seed", &parms->seed, U64, 1, "285920582482");	
   object_get(obj, "randomizeSeed", &randomize, INT, 1, "0");

   if (randomize != 0)
      parms->seed = generateRandomSeed();
   
   return parms;
}

void replicate(TRANSFORM* transform)
{
   timestamp("Start Writing Replicated Restart");

   REPLICATE_PARMS* parms = transform->parms;
   SIMULATE* simulate = transform->parent;
   STATE* state = simulate->system->collection->state;
   unsigned nlocal = simulate->system->nlocal;
   gid_type nglobal = simulate->system->nglobal;

   double* rx = state->rx;
   double* ry = state->ry;
   double* rz = state->rz;
   double* vx = state->vx;
   double* vy = state->vy;
   double* vz = state->vz;

   int nx=parms->nx;
   int ny=parms->ny;
   int nz=parms->nz;
   unsigned nCopies = nx * ny * nz;
   
   char tmp[32];
   OBJECT* obj = (OBJECT*) transform;
   sprintf(tmp, "%f", simulate->time);
   object_get(obj, "time", &transform->time, DOUBLE, 1, tmp);
   transform->time = units_convert(transform->time,"t",NULL);
   sprintf(tmp, "%"PRId64, simulate->loop);
   object_get(obj, "loop", &transform->loop, U64,    1, tmp);
   strcpy(tmp,"transform.");
   sprintf(tmp+strlen(tmp), loopFormat(), transform->loop);
   object_get(obj, "writedir", &transform->writedir, STRING, 1, tmp);
   
   simulate->time = transform->time;
   simulate->loop = transform->loop;
   CreateSnapshotdir(simulate, transform->writedir);
   char filename[512];
   sprintf(filename, "%s/atoms", simulate->snapshotdir);
   PFILE* file = Popen(filename, "w", COMM_LOCAL);
   
   void (*blockWriteFunction) (SIMULATE*simulate, PFILE *file);
   if (simulate->checkpointmode == BINARY_CHECKPOINT)
      blockWriteFunction = collection_writeBLOCK_binary;
   else
      blockWriteFunction = collection_writeBLOCK;

   gid_type* label_new = ddcMalloc(nlocal*sizeof(gid_type));
   double* vx_orig = ddcMalloc(nlocal*sizeof(double));
   double* vy_orig = ddcMalloc(nlocal*sizeof(double));
   double* vz_orig = ddcMalloc(nlocal*sizeof(double));
   gid_type stride = parms->stride;

   gid_type labelHigh[nlocal]; 
   gid_type labelLow[nlocal]; 
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      labelHigh[ii]= state->label[ii]/stride ;
      labelLow[ii] = state->label[ii]%stride ;
      label_new[ii] = nCopies*labelHigh[ii];
      vx_orig[ii] = vx[ii];
      vy_orig[ii] = vy[ii];
      vz_orig[ii] = vz[ii];
   }
   
   gid_type* shift = ddcMalloc(nCopies*sizeof(gid_type));
   for (unsigned ii=0; ii<nCopies; ++ii)
      shift[ii] = ii;
   shuffle(shift, nCopies);

   THREE_MATRIX h = box_get_h(NULL);
   double xsize = h.xx;
   double ysize = h.yy;
   double zsize = h.zz;

   if ( h.xy!=0 || h.xz!=0 || h.yx!=0 || h.yz!=0 || h.zx!=0 || h.zy!=0)
   {
      if (getRank(0) ==0)
	 printf("Replicate transform supports only orthorhombic boxes.\n");

      exit(3);
   }
   h.xx *= nx;
   h.yy *= ny;
   h.zz *= nz;

   box_put(NULL, HO, &h);
   box_put(NULL, HFAC, (THREE_MATRIX *)&I_3x3);
   
   arrayShift(rx, nlocal, (1-nx)*xsize/2.0); 
   arrayShift(ry, nlocal, (1-ny)*ysize/2.0); 
   arrayShift(rz, nlocal, (1-nz)*zsize/2.0); 
   
   collection_writeMode(WREPLICATE);
   collection_nCopies(nCopies);
   for (int ii=0; ii<parms->nx; ++ii)
   {
   for (int jj=0; jj<parms->ny; ++jj)
   {
	 for (int kk=0; kk<parms->nz; ++kk)
	 {
	    if (ii+jj+kk > 0) collection_writeMode(NOWHEADER);
	    for (unsigned ll=0; ll<nlocal; ++ll)
	    {
//	       int iShift = (ii*ny*nz+jj*nz+kk) % nCopies; 
	       int iShift = ((label_new[ll]/nCopies)+ii*ny*nz+jj*nz+kk) % nCopies; 
          gid_type highNew  = (label_new[ll] + shift[iShift]);
	       state->label[ll] = labelLow[ll] + highNew*stride; 
          /*
          if (ll ==-13) 
          {
             printf("label=%16"PRIu64" %16"PRIu64,labelHigh[ll],labelLow[ll]); 
             printf(" %16"PRIu64 "%16"PRIu64"\n",state->label[ll],highNew); 
          }
          */
	       PRAND48_STATE handle =
		  prand48_init(state->label[ll], parms->seed, 0x4321edcb9876llu);
	       vx[ll] = vx_orig[ll] + parms->gamma*gasdev0(&handle);
	       vy[ll] = vy_orig[ll] + parms->gamma*gasdev0(&handle);
	       vz[ll] = vz_orig[ll] + parms->gamma*gasdev0(&handle);
	       if (state->group[ll]->parse != NULL)
		  state->group[ll]->parse(state->group[ll], NULL, ll);
	    }
	    blockWriteFunction(simulate, file);
	    arrayShift(rz, nlocal, zsize);
	 }
	 arrayShift(rz, nlocal, -zsize*nz);
	 arrayShift(ry, nlocal, ysize);
      }
      arrayShift(ry, nlocal, -ysize*ny);
      arrayShift(rx, nlocal, xsize);
   }
   
   Pclose(file);
   if (getRank(0) == 0)
   {
      sprintf(filename, "%s/restart", simulate->snapshotdir);
      FILE* file = fopen(filename, "w");
      fprintf(file, "%s SIMULATE { run_id=0x%08x; loop=%"PRId64"; time=%f ;}\n",
	      simulate->name, simulate->run_id, simulate->loop, simulate->time);
      
      box_write(simulate->system->box, file);
      fprintf(file, "%s COLLECTION { size=%"PRIu64"; files=%s/atoms#;}\n",
	      simulate->system->collection->name,
	      nglobal*nCopies, simulate->snapshotdir);
      fclose(file);
   }
   timestamp("Finish Writing Restart");
   collection_writeMode(WHEADER);
   ddcFree(label_new);
   ddcFree(shift);
}

void arrayShift(double* x, unsigned n, double dx)
{
   for (unsigned ii=0; ii<n; ++ii)
      x[ii] += dx;
}

/** This is a pretty pointless shuffle routine.  In this file it is used
 *  only once to scramble an array with n elements.  If you're going to copy
 *  this to a place where you need real random shuffles be sure to get a
 *  better seed/generator scheme.  In this context the most important
 *  thing is that all tasks have the same shuffle.  Otherwise the
 *  results won't be reproducible on different numbers of tasks. */
void shuffle(gid_type* a, unsigned n)
{
   unsigned short seedVec[3] = {0xf4a2, 0xb58d, 0x90ba};
   for (unsigned ii=0; ii<n-1; ++ii)
   {
      unsigned jj = ii + (unsigned) (erand48(seedVec) * (n - ii));
      gid_type tmp = a[ii];
      a[ii] = a[jj];
      a[jj] = tmp;
   }
}
