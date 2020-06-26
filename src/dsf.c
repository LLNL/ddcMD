#include <stdio.h>
#include <complex.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>

#include "analysis.h"
#include "realsizec.h"
#include "ddcMalloc.h"
#include "simulate.h"
#include "object.h"
#include "species.h"

static void addKvectors(int m, unsigned* nk, THREE_INT** kvec);
int getRank(int);

typedef struct dsf_parms_st
{
   unsigned   nk;
   int        mMax;
   unsigned   nBuf;
   unsigned   nBufMax;
   unsigned*  loop;
   double*    time;
   double complex*   buffer;
   FILE*      file;
   THREE_INT* kvec;
   SPECIES*   species;
} DSF_PARMS;



DSF_PARMS* dsf_parms(ANALYSIS *analysis)
{
   DSF_PARMS* parms = (DSF_PARMS*) ddcMalloc(sizeof(DSF_PARMS)); 
	
   OBJECT* obj = (OBJECT*) analysis;

   char* filename = NULL;
   char* speciesName = NULL;
   int* m = NULL;
   char fileDefault[256];
   int nm = object_getv(obj, "m", (void*)&m, INT,ABORT_IF_NOT_FOUND);
   object_get(obj, "species", &speciesName, STRING, 1, "");
   parms->species = NULL;
   sprintf(fileDefault, "rho_k");
   if (speciesName != NULL)
   {
      parms->species = species_find(NULL, speciesName);
      assert(parms->species != NULL);
      strcat(fileDefault, "_");
      strcat(fileDefault, parms->species->name);
   }
   strcat(fileDefault, ".data");
   object_get(obj, "filename", &filename, STRING, 1, fileDefault);
   
   parms->nk = 0;
   parms->kvec = NULL;
   parms->mMax = 0;
   for (int ii=0; ii<nm; ++ii)
   {
      addKvectors(m[ii], &parms->nk, &parms->kvec);
      if (m[ii] > parms->mMax) parms->mMax = m[ii];
   }
   
   parms->nBuf = 0;
   parms->nBufMax = analysis->outputrate/analysis->eval_rate + 1;
   
   parms->loop   = (unsigned*) ddcMalloc(parms->nBufMax*sizeof(unsigned));
   parms->time   = (double*)   ddcMalloc(parms->nBufMax*sizeof(double));
   parms->buffer = (double complex*)  ddcMalloc(parms->nBufMax*parms->nk*sizeof(double complex));

   if (getRank(0) == 0)
   {
      assert(filename != NULL);
      parms->file = fopen(filename, "a");
      assert(parms->file != NULL);
      
      fprintf(parms->file, "%-8s %16s", "#loop", "time");
      for (unsigned ii=0; ii<parms->nk; ++ii)
      {
         char tmp[256];
         sprintf(tmp, "    (%d,%d,%d)",
                 parms->kvec[ii].x, parms->kvec[ii].y, parms->kvec[ii].z);
         fprintf(parms->file, "%-30s", tmp);
      }
      fprintf(parms->file, "\n");
      fflush(parms->file);
   }
   
   ddcFree(filename);
   ddcFree(speciesName);
   ddcFree(m);

   return parms;
}

void dsf_output(ANALYSIS* analysis)
{
   DSF_PARMS* parms = (DSF_PARMS*) analysis->parms;

   unsigned  nk = parms->nk;
   double complex*  buffer = parms->buffer;
   unsigned* loop = parms->loop;
   double*   time = parms->time;
   FILE*     file = parms->file;
   
   if (getRank(0) == 0)
   {
      for (unsigned ii=0; ii<parms->nBuf; ++ii)
      {
         fprintf(file, "%8.8d %16.6f", loop[ii], time[ii]);
         for (unsigned jj=0; jj<nk; ++jj)
         {
            double rp = creal(buffer[ii*nk+jj]);
            double ip = cimag(buffer[ii*nk+jj]);
            fprintf(file, "   %13.6e %13.6e", rp, ip);
         }
         fprintf(file, "\n");
      }
      fflush(file);
   }
   parms->nBuf = 0;
}


void dsf_eval(ANALYSIS* analysis)
{
   DSF_PARMS* parms = (DSF_PARMS*) analysis->parms;
   unsigned nBuf = parms->nBuf;
   unsigned nBufMax = parms->nBufMax;
   unsigned mMax = parms->mMax;
   unsigned nk = parms->nk;
   THREE_INT* kvec = parms->kvec;
   double complex* buffer = parms->buffer;
   unsigned* loop = parms->loop;
   double*   time = parms->time;
   
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   SYSTEM* sys=simulate->system;
   unsigned nlocal = sys->nlocal;
   double* rx = sys->collection->state->rx;
   double* ry = sys->collection->state->ry;
   double* rz = sys->collection->state->rz;
   double* qi = sys->collection->state->q;
   SPECIES** species = sys->collection->state->species;
   
   if (nBuf >= nBufMax)
      dsf_output(analysis);

   THREE_VECTOR kbasis[3]; 
   box_get(sys->box, RECIP_LATTICEVECTORS, kbasis); 

   // first build up the "basis"-vectors
   static double complex *px=NULL,*py=NULL,*pz=NULL;
   px = (double complex *)ddcRealloc(px, nlocal*(mMax+1)*sizeof(*px)); 
   py = (double complex *)ddcRealloc(py, nlocal*(mMax+1)*sizeof(*py));
   pz = (double complex *)ddcRealloc(pz, nlocal*(mMax+1)*sizeof(*pz));

   gid_type count=0;
   for (unsigned ii=0; ii<nlocal; ++ii)
   {
      if (parms->species && species[ii] != parms->species)
         continue;
      ++count;
      THREE_VECTOR rv; 
      double complex* pxi = px+ii; 
      double complex* pyi = py+ii; 
      double complex* pzi = pz+ii; 
      rv.x=rx[ii]; rv.y=ry[ii]; rv.z=rz[ii]; 
      double complex cx = cexp(I*DOT(kbasis[0], rv));
      double complex cy = cexp(I*DOT(kbasis[1], rv));
      double complex cz = cexp(I*DOT(kbasis[2], rv));
      pxi[0] =  pyi[0] = pzi[0] =  1.0;  
      double complex sx,sy,sz; 
      sx = sy = sz = 1.0; 
      for (unsigned m=1; m<=mMax; m++) 
      {
         pxi[m*nlocal] = (sx *= cx);
         pyi[m*nlocal] = (sy *= cy);
         pzi[m*nlocal] = (sz *= cz);
      }
   }
   
   double complex* localSum = ddcMalloc(nk*sizeof(double complex));
   for (unsigned ii=0; ii<nk; ++ii)
   {
      localSum[ii] = 0;
      double complex* pxi = px+nlocal*abs(kvec[ii].x);
      double complex* pyi = py+nlocal*abs(kvec[ii].y);
      double complex* pzi = pz+nlocal*abs(kvec[ii].z);
      
      for (unsigned jj=0; jj<nlocal; ++jj)
      {
         if (parms->species && species[jj] != parms->species)
            continue;
         double complex xx = pxi[jj];
         double complex yy = pyi[jj];
         double complex zz = pzi[jj];
         if (kvec[ii].x < 0) xx = conj(xx);
         if (kvec[ii].y < 0) yy = conj(yy);
         if (kvec[ii].z < 0) zz = conj(zz);
         localSum[ii] += qi[jj] * xx * yy * zz;
      }
   }

   gid_type countSum=0;
   MPI_Reduce(localSum, buffer+nBuf*nk, 2*nk, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
   MPI_Reduce(&count, &countSum, 1, MPI_GID_TYPE, MPI_SUM, 0, COMM_LOCAL);
   if (countSum > 0)
      for (unsigned ii=0; ii<nk; ++ii)
         buffer[nBuf*nk+ii] /= countSum;
   loop[nBuf] = simulate->loop;
   time[nBuf] = simulate->time;
   ++parms->nBuf;
   ddcFree(localSum);
}

void dsf_close(ANALYSIS* analysis)
{
   DSF_PARMS* parms = (DSF_PARMS*) analysis->parms;
   if (getRank(0) == 0)
   {
      dsf_output(analysis);
      fclose(parms->file);
   }
   
   ddcFree(parms->loop);
   ddcFree(parms->time);
   ddcFree(parms->buffer);
}


void addKvectors(int m, unsigned* nk, THREE_INT** kvec)
{
   unsigned nkMax = *nk;
   *kvec = (THREE_INT*) ddcRealloc(*kvec, nkMax*sizeof(THREE_INT));

   unsigned msq = m*m;
   for (int ii=-m; ii<=m; ++ii)
   {
      unsigned isq = ii*ii;
      for (int jj=-m; jj<=m; ++jj)
      {
         unsigned jsq = jj*jj;
         for (int kk=-m; kk<=m; ++kk)
         {
            unsigned ksq = kk*kk;
            if (kk < 0 ) continue;
            if (kk==0 && jj<0) continue;
            if (kk==0 && jj==0 && ii<=0) continue;

            // testing!!!!!
            if ( ((ii==0) + (jj==0) + (kk==0)) != 2) continue; 
	    
            if (isq+jsq+ksq == msq)
            {
               if (nkMax == *nk)
               {
                  nkMax += 50;
                  *kvec = (THREE_INT*) ddcRealloc(*kvec, nkMax*sizeof(THREE_INT));
               }
               VSET(((*kvec)[*nk]), ii, jj, kk);
               ++(*nk);
            }
         }
      }
   }
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
