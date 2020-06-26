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

// maximum length of the miscInfo string that contains the analysis
// parameter infomation that is written in the header of the file.
#define miscInfoSize 256

enum PC_METHOD {PC_GEOM, PC_GRID, PC_NBR_LIST};

typedef struct paircorrelation_st
{
   char *filename; 
   int nsample;
   unsigned nBins;
   unsigned nspecies; 
   enum PC_METHOD method;
   double delta_r,rmin,rmax,*g;
   double logDelta;
   char miscInfo [miscInfoSize];
   double* binLeft;
   double* binRight;
   int (*whichBin) (double r, const struct paircorrelation_st* parms);
}  PAIRCORRELATION_PARMS; 

static void paircorrelation_eval_geom(ANALYSIS* analysis);
static void paircorrelation_eval_grid(ANALYSIS* analysis);
static void paircorrelation_eval_nbrList(ANALYSIS* analysis);

static int comboIndex(int i, int j, int nspecies);
static void comboReverseIndex(int l, int nspecies, int* iout, int* jout);

static int linearBins(double r, const PAIRCORRELATION_PARMS* parms)
{
   return (r-parms->rmin)/parms->delta_r;
}

static int logBins(double r, const PAIRCORRELATION_PARMS* parms)
{
   return (log10(r)-log10(parms->rmin))/parms->logDelta;
}


PAIRCORRELATION_PARMS *paircorrelation_parms(ANALYSIS *analysis)
{
   analysis->parms = (PAIRCORRELATION_PARMS *) ddcMalloc(sizeof(PAIRCORRELATION_PARMS));
   PAIRCORRELATION_PARMS* parms = analysis->parms;
   OBJECT* obj = (OBJECT*) analysis;

   object_get(obj, "filename", &parms->filename, STRING, 1, "paircorrelation.dat");
   object_get(obj, "length", &parms->nBins, INT, 1, "1");
   object_get(obj, "delta_r", &parms->delta_r, WITH_UNITS, 1, "1","l",NULL);
   object_get(obj, "rmin",    &parms->rmin,    WITH_UNITS, 1, "0","l",NULL);

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
      printf("ERROR: Unrecognized method \"%s\" in pair correlation ANALYSIS object.\n", tmp);
      abortAll(-1);
   }
   
   parms->rmax = parms->rmin + parms->nBins*parms->delta_r; 

   parms->binLeft  = ddcMalloc(parms->nBins * sizeof(double));
   parms->binRight = ddcMalloc(parms->nBins * sizeof(double));
   object_get(obj, "rscale", &tmp, STRING, 1, "normal");
   parms->whichBin = NULL;
   parms->logDelta = 0;
   if (strcasecmp(tmp, "normal") == 0)
   {
      parms->whichBin = linearBins;
      for (unsigned ii=0; ii<parms->nBins; ++ii)
      {
	 parms->binLeft[ii]  = parms->rmin + ii*parms->delta_r;
	 parms->binRight[ii] = parms->rmin + (ii+1)*parms->delta_r;
      }
   }
   else if (strcasecmp(tmp, "log") == 0)
   {
      if (parms->rmin <= 0 && getRank(0) == 0)
      {
	 printf("ERROR: Bad value in pair correlation ANALYSIS object.  When rscale=log\n"
		"rmin must be greater than zero.\n");
	 abortAll(-1);
      }
      parms->whichBin = logBins;
      parms->logDelta = (log10(parms->rmax)-log10(parms->rmin))/(parms->nBins*1.0);
      for (unsigned ii=0; ii<parms->nBins; ++ii)
	 parms->binLeft[ii] = pow(10, log10(parms->rmin) + ii*parms->logDelta);
      for (unsigned ii=0; ii<parms->nBins-1; ++ii)
	 parms->binRight[ii] = parms->binLeft[ii+1];
      parms->binRight[parms->nBins-1] = parms->rmax;
   }
   else if (getRank(0) == 0)
   {
      printf("ERROR: Unrecognized value for rscale keyword \"%s\" in pair correlation ANALYSIS object.\n", tmp);
      abortAll(-1);
   }
   
   species_get(NULL, NSPECIES, (void*) &parms->nspecies);
   parms->g = ddcMalloc(sizeof(double)*parms->nBins*((parms->nspecies*(parms->nspecies+1))/2)); 
   for (unsigned k=0; k<parms->nBins*parms->nspecies*(parms->nspecies+1)/2; k++)
      parms->g[k] = 0.0; 
   parms->nsample = 0; 
   snprintf(parms->miscInfo, miscInfoSize,
	    "rmin = %f Ang; delta_r = %f Ang; length = %d; eval_rate = %d; outputrate = %d;",
	    units_convert(parms->rmin,   NULL, "Angstrom"),
	    units_convert(parms->delta_r, NULL, "Angstrom"),
	    parms->nBins, analysis->eval_rate, analysis->outputrate);

   return parms; 
}

void paircorrelation_clear(ANALYSIS* analysis)
{
   PAIRCORRELATION_PARMS* parms = (PAIRCORRELATION_PARMS *)analysis->parms;
   for (unsigned k=0; k<parms->nBins*parms->nspecies*(parms->nspecies+1)/2; k++)
      parms->g[k] = 0.0; 
   parms->nsample=0; 
}


void paircorrelation_eval(ANALYSIS* analysis)
{
   PAIRCORRELATION_PARMS* parms = (PAIRCORRELATION_PARMS *)analysis->parms;
   switch (parms->method)
   {
     case PC_NBR_LIST:
      paircorrelation_eval_nbrList(analysis);
      break;
     case PC_GEOM:
      paircorrelation_eval_geom(analysis);
      break;
     case PC_GRID:
      paircorrelation_eval_grid(analysis);
      break;
     default:
      assert(1==0);
   }
}

void paircorrelation_eval_nbrList(ANALYSIS *analysis)
{
   SIMULATE *simulate; 
   SYSTEM *sys; 
   STATE *state; 
   PAIRCORRELATION_PARMS *parms ; 
   double r2,rmax2, minspan,R2cut;
   double *rx,*ry,*rz,x0,y0,z0,x,y,z; 
   int i,j,k,Ti,Tj,nlocal;
   PAIRS *pij; 
   int *type; 
   parms = (PAIRCORRELATION_PARMS *)analysis->parms; 
   simulate =(SIMULATE *)analysis->parent; 
   sys = simulate->system; 
   state = sys->collection->state; 
   nlocal = sys->nlocal; 
   rx = state->rx; 
   ry = state->ry; 
   rz = state->rz; 
   type = state->atomtype;
   box_get(NULL,MINSPAN,&minspan);
   R2cut = 0.25*minspan*minspan;
   rmax2 = parms->rmax*parms->rmax;
   parms->nsample +=1; 
   NBR *nbr = sys->neighbor;
   pij=nbr->pairs; 

   unsigned ns = parms->nspecies; 
   unsigned np = (ns*(ns+1))/2;
   unsigned nbl = parms->nBins*np;
   double nBonds[nbl];
   double nAtoms[ns];
   for (unsigned ii=0; ii<nbl; ++ii) nBonds[ii] = 0;
   for (unsigned ii=0; ii<ns; ++ii) nAtoms[ii] = 0;
   for (i = 0; i < nlocal; i++)
   {
      NPARTICLE *pi = nbr->particles + i;
      x0 = rx[i]; y0 = ry[i]; z0 = rz[i]; 
      Ti = type[i]&0xffff;
      nAtoms[Ti]++;
      for (k = 0;k<pi->nb;k++)
      {
	 j=pij->j ;
	 Tj = type[j]&0xffff;
	 x = x0 - rx[j]; y = y0 - ry[j]; z = z0 - rz[j]; r2 = x*x+y*y+z*z;
	 if (r2 > R2cut) {nearestImage_fast(&x, &y, &z); r2 = x*x + y*y + z*z;}
	 if (r2 < rmax2)
	 {
	    int k = parms->whichBin(sqrt(r2), parms);
	    int l = comboIndex(Ti, Tj, parms->nspecies);
	    k += parms->nBins*l; 
	    if (Ti == Tj)
	       nBonds[k]+=2.0;
	    else
	       nBonds[k]+=1.0;
	 }
	 pij++; 
      }
   }
   double nAtomsSum[ns];
   MPI_Allreduce(nAtoms, nAtomsSum, ns, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   for (unsigned ii=0; ii<ns; ++ii)
      for (unsigned jj=ii; jj<ns; ++jj)
      {
	 int l = comboIndex(ii, jj, parms->nspecies);
	 double recipNiNj = 1.0/(nAtomsSum[ii]*nAtomsSum[jj]);
	 for (unsigned kk=0; kk<parms->nBins; ++kk)
	    nBonds[kk+parms->nBins*l] *= recipNiNj;
      }
   for (unsigned ii=0; ii<nbl; ++ii)
      parms->g[ii] += nBonds[ii];
   
}

void paircorrelation_eval_grid(ANALYSIS *analysis)
{
   PAIRCORRELATION_PARMS* parms = (PAIRCORRELATION_PARMS* ) analysis->parms; 
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
   double nBonds[nbl];
   double nType[ns];
   for (unsigned ii=0; ii<nbl; ++ii) nBonds[ii] = 0;
   for (unsigned ii=0; ii<ns; ++ii) nType[ii] = 0;

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

	    int k = parms->whichBin(r, parms);
	    int ci = comboIndex(Ti, Tj, parms->nspecies);
	    k += parms->nBins*ci; 
	    if (Ti == Tj)
	       nBonds[k]+=2.0;
	    else
	       nBonds[k]+=1.0;
	 }      
      }
   }

   double nTypeSum[ns];
   MPI_Allreduce(nType, nTypeSum, ns, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   for (unsigned ii=0; ii<ns; ++ii)
      for (unsigned jj=ii; jj<ns; ++jj)
      {
	 int l = comboIndex(ii, jj, parms->nspecies);
	 double recipNiNj = 1.0/(nTypeSum[ii]*nTypeSum[jj]);
	 for (unsigned kk=0; kk<parms->nBins; ++kk)
	 {
	    nBonds[kk+parms->nBins*l] *= recipNiNj;
	 }
      }
   for (unsigned ii=0; ii<nbl; ++ii)
      parms->g[ii] += nBonds[ii];

   pfs_destroy(&pairFinder);
   heapFree(gBlk);
   heapFree(tBlk);
   heapFree(rBlk);
}

void paircorrelation_eval_geom(ANALYSIS *analysis)
{
   PAIRCORRELATION_PARMS* parms = (PAIRCORRELATION_PARMS* ) analysis->parms; 
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
   
   unsigned ns = parms->nspecies; 
   unsigned np = (ns*(ns+1))/2;
   unsigned nbl = parms->nBins*np;
   double nBonds[nbl];
   double nType[ns];
   for (unsigned ii=0; ii<nbl; ++ii) nBonds[ii] = 0;
   for (unsigned ii=0; ii<ns; ++ii) nType[ii] = 0;

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
	    assert(r <= rCut);

	    int k = parms->whichBin(r, parms);
	    int ci = comboIndex(Ti, Tj, parms->nspecies);
	    k += parms->nBins*ci; 
	    if (Ti == Tj)
	       nBonds[k]+=2.0;
	    else
	       nBonds[k]+=1.0;
	 }      
      }
   }

   double nTypeSum[ns];
   MPI_Allreduce(nType, nTypeSum, ns, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   for (unsigned ii=0; ii<ns; ++ii)
      for (unsigned jj=ii; jj<ns; ++jj)
      {
	 int l = comboIndex(ii, jj, parms->nspecies);
	 double recipNiNj = 1.0/(nTypeSum[ii]*nTypeSum[jj]);
	 for (unsigned kk=0; kk<parms->nBins; ++kk)
	 {
	    nBonds[kk+parms->nBins*l] *= recipNiNj;
	 }
      }
   for (unsigned ii=0; ii<nbl; ++ii)
      parms->g[ii] += nBonds[ii];

   pfg_destroy(&pairFinder);
   heapFree(gBlk);
   heapFree(tBlk);
   heapFree(rBlk);
}


void  paircorrelation_output(ANALYSIS *analysis)
{
   PAIRCORRELATION_PARMS* parms = (PAIRCORRELATION_PARMS *)analysis->parms; 
   SIMULATE *simulate =(SIMULATE *)(analysis->parent); 
   int ns = parms->nspecies; 
   int np = (ns*(ns+1))/2;
   if (parms->nsample == 0 )
   {
      paircorrelation_clear(analysis);
      return ;
   }
   double* g = ddcMalloc(sizeof(double)*np*parms->nBins); 
   MPI_Reduce(parms->g, g, np*parms->nBins, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
   
   CreateSnapshotdir(simulate, NULL);
   if (getRank(0) == 0 ) 
   {
      char* sName[ns];
      for (int ii=0; ii<ns; ++ii)
	 sName[ii] = species_by_index(NULL, ii)->name;
      char filename[1024];
      snprintf(filename, 1024,"%s/%s", simulate->snapshotdir,parms->filename);
      FILE* file = fopen(filename, "w");
      double volume=0;
      box_get(NULL,VOLUME,&volume);
      double s = volume/(parms->nsample);
      for (unsigned  k=0;k<parms->nBins;k++)
      {
	 double r0 =  parms->binLeft[k];
	 double r1 =  parms->binRight[k];
	 double dv = 4.0*M_PI /3.0 * (r1*r1*r1-r0*r0*r0);
	 for (int l=0; l<np; l++)
	    g[k+parms->nBins*l] *= s/dv;
      }
      
      fprintf(file, "# %s\n", parms->miscInfo);
      fprintf(file, "# nsample = %d;\n", parms->nsample);
      fprintf(file, "# r(Ang) ");
      for (int jj=0; jj<np; ++jj)
      {
	 int Ti=-1, Tj=-1;
	 comboReverseIndex(jj, parms->nspecies, &Ti, &Tj);
	 fprintf(file, "%s-%s ", sName[Ti], sName[Tj]);
      }
      fprintf(file, "\n");
      
      for (unsigned k=0;k<parms->nBins;k++)  
      {
	 double binCenter = 0.5*(parms->binLeft[k] + parms->binRight[k]);
	 fprintf(file,"%f ", units_convert(binCenter, NULL, "Angstrom"));
	 for (int l=0;l<np;l++)
	    fprintf(file,"%e ",g[k+parms->nBins*l]); 
	 fprintf(file,"\n");
      }
      fclose(file);
   }
   paircorrelation_clear(analysis);
   ddcFree(g);
}

int comboIndex(int i, int j, int nspecies)
{
   int max_ij = MAX(i, j); 
   int min_ij = MIN(i, j); 
   return ((max_ij-min_ij) + nspecies*min_ij-(min_ij*(min_ij-1))/2);
}

/** Yes, this is dreadfully inefficient.  But, this implementation is
 * guaranteed to stay in sync with any possible changes to the
 * comboIndex function.  It only happens once per run and the number of
 * species is practically certain to be small.  Hence the search space
 * is actually quite small and this will not take any noticable amount
 * of time.  */
void comboReverseIndex(int l, int nspecies, int* iout, int* jout)
{
   for (int ii=0; ii<nspecies; ++ii)
      for (int jj=ii; jj<nspecies; ++jj)
      {
	 if (comboIndex(ii, jj, nspecies) == l)
	 {
	    *iout = ii;
	    *jout = jj;
	    return;
	 }
      }
   printf("l=%d\n",l);
   assert(0==1);
}

