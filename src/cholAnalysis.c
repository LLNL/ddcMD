
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>

#include "analysis.h"
#include "realsizec.h"
#include "ddcMalloc.h"
#include "simulate.h"
#include "object.h"
#include "error.h"
#include "three_algebra.h"
#include "expandbuffer.h"
#include "io.h"
#include "units.h"
#include "codata.h"
#include "preduce.h"
#include "bioCharmm.h"
#include "bioCharmmParms.h"

THREE_VECTOR mBond(THREE_VECTOR r0, THREE_VECTOR r1);
typedef struct cholAnalysis_parms_st
{
  char *potentialName;
  FILE *datafile; 
  CHARMMPOT_PARMS *potParms; 
  char *filename; 
  double rmin,rmax,delta; 
  int nbins; 
  int *cnt;
  double dR1Min,dR1Max,dR1Ave; 
  double dR5Min,dR5Max,dR5Ave; 
} CHOLANALYSIS_PARMS;

int getRank(int);

CHOLANALYSIS_PARMS* cholAnalysis_parms(ANALYSIS *analysis)
{
   OBJECT* obj = (OBJECT*) analysis;
   char *dataFilename; 
   CHOLANALYSIS_PARMS* parms = (CHOLANALYSIS_PARMS*) ddcMalloc(sizeof(CHOLANALYSIS_PARMS)); 
	object_get(obj, "potentialName", &parms->potentialName, STRING, 1, "martini");
	object_get(obj, "filename", &parms->filename, STRING, 1, "cholAnalysis.distn");
	object_get(obj, "dataFilename", &dataFilename, STRING, 1, "cholAnalysis.data");
   parms->datafile = fopen(dataFilename,"a"); 
   POTENTIAL *potential = (POTENTIAL*) object_longFind(parms->potentialName,"POTENTIAL");
   assert(potential != NULL); 
   assert(strcmp(potential->type ,"MARTINI") ==0 ); 
   parms->potParms=(CHARMMPOT_PARMS *)potential->parms; 
   assert(parms->potParms != NULL); 
   object_get(obj, "delta", &parms->delta, WITH_UNITS, 1, "0.1","l",NULL);
   object_get(obj, "rmin",    &parms->rmin,    WITH_UNITS, 1, "0","l",NULL);
   object_get(obj, "rmax",    &parms->rmax,    WITH_UNITS, 1, "0","l",NULL);
   parms->nbins = lrint((parms->rmax-parms->rmin)/parms->delta); 
   parms->delta = (parms->rmax-parms->rmin)/parms->nbins; 
   parms->cnt = (int *)ddcMalloc(sizeof(int)*2*parms->nbins); 
   for (int i=0;i<2*parms->nbins;i++) parms->cnt[i]=0; 
   parms->dR1Min = parms->dR5Min = 1e+300; 
   parms->dR1Max = parms->dR5Max = -1e+300; 
   parms->dR1Ave = parms->dR5Ave = 0; 
	
	

   //object_get(obj, "kmax", &parms->kmax,DOUBLE,1,"1.0");
   
   return parms;
}

void cholAnalysis_output(ANALYSIS* analysis)
{
   CHOLANALYSIS_PARMS* parms = (CHOLANALYSIS_PARMS*) analysis->parms;
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   double lc = units_convert(1.0, NULL, "Angstrom");
	CreateSnapshotdir(simulate, NULL);
   FILE *file=NULL; 
   int cnt1=0; 
   int cnt3=0; 
   for (int i=0;i<parms->nbins;i++)  
   {
      cnt1 += parms->cnt[i];
      cnt3 += parms->cnt[i+parms->nbins];
   }
   printf("# cnt1=%d cnt3=%d\n",cnt1,cnt3); 
	if (getRank(0) == 0 ) 
	{
		char filename[1024]; 
		snprintf(filename, 1023,"%s/%s", simulate->snapshotdir,parms->filename);
		file = fopen(filename, "w");
      parms->dR1Ave /= cnt1; 
      parms->dR5Ave /= cnt1; 
      fprintf(parms->datafile,"%ld %f %f %f %f %f %f %f\n",simulate->loop,simulate->time,parms->dR1Min*lc,parms->dR1Max*lc,parms->dR1Ave*lc,parms->dR5Min*lc,parms->dR5Max*lc,parms->dR5Ave*lc);
      fflush(parms->datafile); 
   }
   parms->dR1Min = parms->dR5Min = 1e+300; 
   parms->dR1Max = parms->dR5Max = -1e+300; 
   parms->dR1Ave = parms->dR5Ave = 0; 
   for (int i=0;i<parms->nbins;i++)  
   {
      double r = (parms->rmin+(i+0.5)*parms->delta);
      fprintf(file," %e %e %e\n",r*lc,parms->cnt[i]/lc*(1.0/(cnt1*parms->delta)), parms->cnt[i+parms->nbins]/lc*(1.0/(cnt3*parms->delta))); 
   }
   fclose(file); 
   for (int i=0;i<2*parms->nbins;i++) parms->cnt[i]=0; 
}
void cholAnalysis_eval(ANALYSIS* analysis)
{
   CHOLANALYSIS_PARMS* parms = (CHOLANALYSIS_PARMS*) analysis->parms;
   
   SIMULATE* simulate =(SIMULATE *)analysis->parent; 
   SYSTEM* sys=simulate->system;
   //unsigned nlocal = sys->nlocal;
   SETLIST residueSet = parms->potParms->residueSet; 
   GID_ORDER *gidOrder=parms->potParms->gidOrder;    
   LISTNODE *residueList=residueSet.list; 
   double* rx = sys->collection->state->rx;
   double* ry = sys->collection->state->ry;
   double* rz = sys->collection->state->rz;
   //SPECIES**species = sys->collection->state->species;
   //double lc = units_convert(1.0, NULL, "Angstrom");
   int index =0; 
    for (int i = 0; i < residueSet.listSize; i++)
    {
       RESI_CONN *resiConn = residueList[i].resiConn;
       char* name = resiConn->resName;
       if (strcmp(name,"CHOL") == 0) 
       { 
       THREE_VECTOR r[resiConn->atomListSize]; 
       for (int j =0; j<resiConn->atomListSize; j++)
       {
         int k = gidOrder[index+j].id;
         VSET(r[j],rx[k],ry[k],rz[k]);
       }
       THREE_VECTOR A = mBond(r[0],r[1]); 
       THREE_VECTOR B = mBond(r[0],r[2]); 
       THREE_VECTOR C = mBond(r[0],r[3]); 
       THREE_VECTOR D = mBond(r[4],r[5]); 
       THREE_VECTOR E = mBond(r[4],r[3]); 
       THREE_VECTOR F = mBond(r[4],r[6]); 
       //THREE_VECTOR x1,x2,x3; 
       THREE_VECTOR x1,x3; 
       CROSS(x1,B,C); 
       double dR1 = DOT(x1,A)/sqrt(VSQ(x1)); 
       if (dR1 < parms->dR1Min) parms->dR1Min=dR1; 
       if (dR1 > parms->dR1Max) parms->dR1Max=dR1; 
       parms->dR1Ave += dR1; 
       CROSS(x3,E,F); 
       double dR5 = -DOT(x3,D)/sqrt(VSQ(x3)); 
       if (dR5 < parms->dR5Min) parms->dR5Min=dR5; 
       if (dR5 > parms->dR5Max) parms->dR5Max=dR5; 
       parms->dR5Ave += dR5; 
//       printf("delta %f %f\n",dR1*lc,dR5*lc); 
       int ibin1 =  MIN(MAX((dR1-parms->rmin)/parms->delta,0),parms->nbins-1); 
       int ibin3 =  MIN(MAX((dR5-parms->rmin)/parms->delta,0),parms->nbins-1); 
       parms->cnt[ibin1]++;
       parms->cnt[ibin3+parms->nbins]++;
       }
       index += resiConn->atomListSize;
    }
}
void cholAnalysis_close(ANALYSIS* analysis)
{
}
THREE_VECTOR mBond(THREE_VECTOR r0, THREE_VECTOR r1) 
{
   THREE_VECTOR r; 
   VOP2(r,=,r1,-,r0); 
   nearestImage(&r.x, &r.y, &r.z);
   return r; 
}
void lookat(int loop,unsigned gid, int location, THREE_VECTOR *r)
{
   THREE_VECTOR A = mBond(r[0],r[1]); 
   THREE_VECTOR B = mBond(r[0],r[2]); 
   THREE_VECTOR C = mBond(r[0],r[3]); 
   THREE_VECTOR D = mBond(r[6],r[5]); 
   THREE_VECTOR E = mBond(r[6],r[4]); 
   THREE_VECTOR F = mBond(r[6],r[3]); 
   THREE_VECTOR x; 
   CROSS(x,B,C); 
   double d0 = DOT(x,A)/sqrt(VSQ(x)); 
   CROSS(x,E,F); 
   double d1 = DOT(x,D)/sqrt(VSQ(x)); 
   printf("delta %d %d %d %f %f\n",gid,loop,location,d0,d1); 

}
/* Local Variables: */
/* tab-width: 3 */
/* End: */
