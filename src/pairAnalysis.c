#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "three_algebra.h"
#include "object.h"
#include "pio.h"
#include "ddc.h"
#include "utilities.h"
#include "crc32.h"
#include "simulate.h"
#include "system.h"
#include "expandbuffer.h"
#include "ddcMalloc.h"
#include "io.h"
#include "box.h"
#include "units.h"
#include "neighbor.h"
#include "heap.h"
#include "dataExchange.h"
#include "dataExchangeUtils.h"
#include "pairIterator.h"
#include "pairFinder.h"
#include "pairFinderGeom.h"
#include "mpiUtils.h"
#include "preduce.h"
#include "gid.h"

// maximum length of the miscInfo string that contains the analysis
// parameter infomation that is written in the header of the file.
#define miscInfoSize 256

enum PC_METHOD {PC_GEOM, PC_GRID, PC_NBR_LIST};

typedef struct pair_st
{
  gid_type gid[2];
} PAIRTYPE;
typedef struct pairAnalysis_st
{
   char *filename; 
   unsigned nspecies; 
   double rmax; 
   enum PC_METHOD method;
   PAIRTYPE *pairs;
}  PAIRANALYSIS_PARMS; 

static void pairAnalysis_eval_geom(ANALYSIS* analysis);
static void pairAnalysis_eval_grid(ANALYSIS* analysis);
static void pairAnalysis_eval_nbrList(ANALYSIS* analysis);



PAIRANALYSIS_PARMS *pairAnalysis_parms(ANALYSIS *analysis)
{
   analysis->parms = (PAIRANALYSIS_PARMS *) ddcMalloc(sizeof(PAIRANALYSIS_PARMS));
   PAIRANALYSIS_PARMS* parms = analysis->parms;
   OBJECT* obj = (OBJECT*) analysis;

   object_get(obj, "filename", &parms->filename, STRING, 1, "pairAnalysis.dat");
   object_get(obj, "rmax",    &parms->rmax,    WITH_UNITS, 1, "0","l",NULL);

   char* tmp;
   object_get(obj, "method", &tmp, STRING, 1, "geom");
   if (strcasecmp(tmp, "geom") == 0)
      parms->method = PC_GEOM;
   else if (strcasecmp(tmp, "grid") == 0)
      parms->method = PC_GRID;
   else if (strcasecmp(tmp, "neighborList") == 0)
      parms->method = PC_NBR_LIST;
   else if (getRank(0) == 0)
   {
      printf("ERROR: Unrecognized method \"%s\" in pair ANALYSIS object.\n", tmp);
      abortAll(-1);
   }
   
   species_get(NULL, NSPECIES, (void*) &parms->nspecies);
   return parms; 
}

void pairAnalysis_clear(ANALYSIS* analysis)
{
   PAIRANALYSIS_PARMS* parms = (PAIRANALYSIS_PARMS *)analysis->parms;
}


void pairAnalysis_eval(ANALYSIS* analysis)
{
   PAIRANALYSIS_PARMS* parms = (PAIRANALYSIS_PARMS *)analysis->parms;
   switch (parms->method)
   {
     case PC_NBR_LIST:
      pairAnalysis_eval_nbrList(analysis);
      break;
     case PC_GEOM:
      //pairAnalysis_eval_geom(analysis);
      break;
     case PC_GRID:
      //pairAnalysis_eval_grid(analysis);
      break;
     default:
      assert(1==0);
   }
}

void pairAnalysis_eval_nbrList(ANALYSIS *analysis)
{
   SIMULATE *simulate; 
   SYSTEM *sys; 
   PAIRANALYSIS_PARMS *parms ; 
   double rmax2, minspan,R2cut;
   double *rx,*ry,*rz,x0,y0,z0,x,y,z; 
   PAIRS *pij; 
   int *type; 
   parms = (PAIRANALYSIS_PARMS *)analysis->parms; 
   simulate =(SIMULATE *)analysis->parent; 
   sys = simulate->system; 
   STATE *state = sys->collection->state; 
   rx = state->rx; 
   ry = state->ry; 
   rz = state->rz; 
   SPECIES **species = state->species; 
   type = state->atomtype;
   box_get(NULL,MINSPAN,&minspan);
   R2cut = 0.25*minspan*minspan;
   rmax2 = parms->rmax*parms->rmax;
   NBR *nbr = sys->neighbor;
   pij=nbr->pairs; 
   NPARTICLE *particles = sys->neighbor->particles;

   int cnt =0; 
   for (int i = 0; i < sys->nlocal; i++)
   {
      int si = species[i]->index;
      PAIRS *pij = particles[i].ifirst[0];
      //NPARTICLE *pi = nbr->particles + i;
      x0 = rx[i]; y0 = ry[i]; z0 = rz[i]; 
      int flag = 0; 
      //if (pij != NULL) 
      //{
      //   flag = 1; 
      //    printf("%d : ",i); fflush(stdout); 
      //}
      while (pij  != NULL)
      {
         RPAIR *p = pij->p;
         int j=pij->j ;
         int sj = species[j]->index;
         double r2 = 0; 
         x = x0 - rx[j]; y = y0 - ry[j]; z = z0 - rz[j]; r2 = x*x+y*y+z*z;
         if (r2 > R2cut) {nearestImage_fast(&x, &y, &z); r2 = x*x + y*y + z*z;}
         //printf(" %d %e",j,sqrt(r2)); fflush(stdout); 
         if (r2 < rmax2)
         {
          cnt++; 
         }
         pij = pij->ilink;
      }
      //if (flag) {printf("\n");  fflush(stdout); }
   }
   printf("cnt=%d\n",cnt); 
}
#if 0 
 for (unsigned i = 0; i < sys->nlocal; i++)
   {
      int si = species[i]->index;
      PAIRS *pij = particles[i].ifirst[IndexRmax];
      while (pij  != NULL)
      {
         RPAIR *p = pij->p;
         double r = p->r;
         int j=pij->j;
         int sj = species[j]->index;
         int sij = sj+sys->nspecies*si;
         if (r < parms->ppRcut )
         {
            double dvdr=0.0;
            double vij =  parms->pp[sij]->func(parms->pp[sij]->parms,r,&dvdr) ;
            dvdr /= r;
            er += vij;
            //         p->e += vij;
            p->e = 0;
            THREE_VECTOR fp;
            fp.x = -dvdr*p->x;
            fp.y = -dvdr*p->y;
            fp.z = -dvdr*p->z;
            fx[i]+=fp.x ;
            fy[i]+=fp.y ;
            fz[i]+=fp.z ;
            fx[j]-=fp.x ;
            fy[j]-=fp.y ;
            fz[j]-=fp.z ;
         }
         pij = pij->ilink;
      }
   }


void pairAnalysis_eval_grid(ANALYSIS *analysis)
{
   PAIRANALYSIS_PARMS* parms = (PAIRANALYSIS_PARMS* ) analysis->parms; 
   parms->nsample += 1; 
   double rCut = parms->rmax;
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   SYSTEM* sys=simulate->system; 
   unsigned nLocal = sys->nlocal;
   STATE* state = sys->collection->state;

   DATA_EXCHANGE_PARMS dxParms = mkRemoteAtomDep(rCut);
   unsigned nRemote = dep_nRemote(dxParms);
   unsigned nAtoms = nLocal+nRemote;

   unsigned rBlk, tBlk, gBlk;
   THREE_VECTOR* rr = heapGet(&rBlk);
   heapEndBlock(rBlk, nAtoms*sizeof(THREE_VECTOR));
   int* atomType = heapGet(&tBlk);
   heapEndBlock(tBlk, nAtoms*sizeof(int));
   gid_type* gid = heapGet(&gBlk);
   heapEndBlock(gBlk, nAtoms*sizeof(gid_type));

   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      rr[ii].x = state->rx[ii];
      rr[ii].y = state->ry[ii];
      rr[ii].z = state->rz[ii];
      atomType[ii] = state->atomtype[ii];
      gid[ii] = state->label[ii];
   }

   dxParms.width = sizeof(THREE_VECTOR);
   dep_exchange(rr, rr+nLocal, &dxParms);
   dxParms.width = sizeof(int);
   dep_exchange(atomType, atomType+nLocal, &dxParms);
   dxParms.width = sizeof(gid_type);
   dep_exchange(gid, gid+nLocal, &dxParms);

   dep_free(&dxParms);

   //make atoms compact;
   THREE_VECTOR r0 = rr[0];
   for (unsigned ii=0; ii<nAtoms; ++ii)
   {
      rr[ii].x -= r0.x;
      rr[ii].y -= r0.y;
      rr[ii].z -= r0.z;
      nearestImage_fast(&(rr[ii].x), &(rr[ii].y), &(rr[ii].z));
   }

   PAIR_FINDER_SIMPLE pairFinder = pfs_create(rr, nAtoms, nearestImage_fast, rCut);

   unsigned ns = parms->nspecies; 
   unsigned np = (ns*(ns+1))/2;
   unsigned nbl = parms->nBins*np;
   double nType[ns];

   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      unsigned Ti = atomType[ii] & 0xffff;
      ++nType[Ti];
      for (PAIR_ITERATOR iter=pfs_newIter(ii, &pairFinder);
            iter.ii!=iter.jj; pairIter_advance(&iter))
      {
         if (gid[ii] < gid[iter.jj])
         {
            unsigned Tj = atomType[iter.jj] & 0xffff;
            double r = sqrt(iter.r2);
            assert(r <= rCut);

         }      
      }
   }

   pfs_destroy(&pairFinder);
   heapFree(gBlk);
   heapFree(tBlk);
   heapFree(rBlk);
}

void pairAnalysis_eval_geom(ANALYSIS *analysis)
{
   PAIRANALYSIS_PARMS* parms = (PAIRANALYSIS_PARMS* ) analysis->parms; 
   parms->nsample += 1; 
   double rCut = parms->rmax;
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   SYSTEM* sys=simulate->system; 
   unsigned nLocal = sys->nlocal;
   STATE* state = sys->collection->state;

   // Build a communications routing table that will exchange data for
   // all atoms within distance rCut of any local atom.  
   DATA_EXCHANGE_PARMS dxParms = mkRemoteAtomDep(rCut);
   unsigned nRemote = dep_nRemote(dxParms);
   unsigned nAtoms = nLocal+nRemote;

   unsigned rBlk, tBlk, gBlk;
   THREE_VECTOR* rr = heapGet(&rBlk);
   heapEndBlock(rBlk, nAtoms*sizeof(THREE_VECTOR));
   int* atomType = heapGet(&tBlk);
   heapEndBlock(tBlk, nAtoms*sizeof(int));
   gid_type* gid = heapGet(&gBlk);
   heapEndBlock(gBlk, nAtoms*sizeof(gid_type));

   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      rr[ii].x = state->rx[ii];
      rr[ii].y = state->ry[ii];
      rr[ii].z = state->rz[ii];
      atomType[ii] = state->atomtype[ii];
      gid[ii] = state->label[ii];
   }

   dxParms.width = sizeof(THREE_VECTOR);
   dep_exchange(rr, rr+nLocal, &dxParms);
   dxParms.width = sizeof(int);
   dep_exchange(atomType, atomType+nLocal, &dxParms);
   dxParms.width = sizeof(gid_type);
   dep_exchange(gid, gid+nLocal, &dxParms);

   dep_free(&dxParms);

   unsigned domainId = simulate->ddc->domain_id;
   const THREE_VECTOR* domainCenterP = &(simulate->ddc->domains.domains[domainId].center);
   PARTICLESET* ps = NULL;
   ps = ParticleSet(
         ps, &rr[0].x, &rr[0].y, &rr[0].z, sizeof(THREE_VECTOR),
         NULL, 0, NULL, 0, &nAtoms, &sys->nlocal, domainCenterP);

   PAIR_FINDER_GEOM pairFinder = pfg_create(ps, nearestImage_fast, rCut,box_getBox(NULL));


   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      unsigned Ti = atomType[ii] & 0xffff;
      ++nType[Ti];
      for (PAIR_ITERATOR iter=pfg_newIter(ii, &pairFinder);
            iter.ii!=iter.jj; pairIter_advance(&iter))
      {
         if (gid[ii] < gid[iter.jj])
         {
            unsigned Tj = atomType[iter.jj] & 0xffff;
            double r = sqrt(iter.r2);
            if (r <= rCut)
            {
            }

         }      
      }
   }

   pfg_destroy(&pairFinder);
   heapFree(gBlk);
   heapFree(tBlk);
   heapFree(rBlk);
}
#endif


void  pairAnalysis_output(ANALYSIS *analysis)
{
   PAIRANALYSIS_PARMS* parms = (PAIRANALYSIS_PARMS *)analysis->parms; 
   SIMULATE *simulate =(SIMULATE *)(analysis->parent); 
   if (0) {
      pairAnalysis_clear(analysis);
      return ;
   }
   CreateSnapshotdir(simulate, NULL);
   if (getRank(0) == 0 ) 
   {
      char filename[1024];
      snprintf(filename, 1024,"%s/%s", simulate->snapshotdir,parms->filename);
      FILE* file = fopen(filename, "w");
      fclose(file);
   }
   pairAnalysis_clear(analysis);
}

