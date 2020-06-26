#include "pair.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <sys/time.h>
#include <mpi.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "utilities.h"
#include "system.h"
#include "ddcMalloc.h"
#include "ddc.h"
#include "codata.h"
#include "units.h"
#include "mpiUtils.h"
#include "binProcess.h"
#include "pairProcessGPU.h"
#include "HAVEGPU.h"
/*
typedef struct pair_parms_str
{
	POTENTIAL *potential; 
	unsigned nspecies; 
	char **species_list; 
	double *rcutoff;
	double rmax;
	double(*fcn)(void *,double, double *); 
	void **parms; 
} PAIR_PARMS;

typedef struct lj_struct 
{ 
	double rcut,eps,sigma;
   double shift;
}  LJ_PARMS;
*/


void *table_parms(POTENTIAL *potential);

double LennardJones(LJ_PARMS *ljparms,double r,double *dvdr) 
{
	double sigma = ljparms->sigma;
	double eps = ljparms->eps; 
	double sigma_r = sigma/r; 
	double s2  = sigma_r*sigma_r; 
	double s4  = s2*s2; 
	double s6  = s4*s2; 
	double s12 = s6*s6; 
	*dvdr = 4.0*eps*(-12.0*s12+6.0*s6)/r; 
	double e = 4.0*eps * (s12 - s6)+ljparms->shift; 
   
	return e; 
}

void LennardJones_setShift(LJ_PARMS *ljparms)
{
   double sigma = ljparms->sigma;
   double eps = ljparms->eps;
   double sigma_r = sigma/ljparms->rcut;
   double s2  = sigma_r*sigma_r;
   double s4  = s2*s2;
   double s6  = s4*s2;
   double s12 = s6*s6;
   ljparms->shift = -4.0*eps * (s12 - s6);
}

LJ_PARMS **LennardJones_parms(PAIR_PARMS *parms)
{
	POTENTIAL *potential=parms->potential; 
	int nspecies = parms->nspecies; 
	double cutoff; 
	LJ_PARMS **ljParms = (LJ_PARMS **)ddcMalloc(sizeof(LJ_PARMS*)*nspecies*nspecies);
   object_get((OBJECT*)potential, "cutoff", &cutoff, DOUBLE, 1, "2.5");
	parms->fcn=(double(*)(void *,double,double *))LennardJones;
	parms->rmax = 0.0; 
	SYSTEM *sys= (SYSTEM*)potential->parent; 
	parms->nspecies=nspecies; 
   
   int list[(nspecies*nspecies+nspecies)/2][2];
   int nlist =0; 
	for (int a = 0;a<nspecies;a++) 
   {  
      list[nlist][0] = a; 
      list[nlist][1] = a; 
      nlist++;
   }
	for (int a = 0;a<nspecies;a++)
	for (int b = 0;b<a;b++)
   {
      
      list[nlist][0]=a;
      list[nlist][1]=b;
      nlist++;
   }
   LJ_PARMS *buffer = ddcMalloc(nlist*sizeof(LJ_PARMS));
   for (int k = 0;k<nlist;k++)
   {
      int a = list[k][0]; 
      int b = list[k][1]; 

      int ab = a+b*nspecies;
      int ba = b+a*nspecies;
      int aa = a+a*nspecies;
      int bb = b+b*nspecies;
      SPECIES *speciesA = sys->species[a]; 
      SPECIES *speciesB = sys->species[b]; 
      int lkey = strlen(speciesA->name)+strlen(speciesB->name);

      char keyword[lkey+16]; 
      double eps,sigma; 
      sprintf(keyword,  "eps_%s-%s",speciesA->name,speciesB->name);
      int flag1=object_testforkeyword((OBJECT*) potential,keyword);
      if (flag1) object_get_with_units((OBJECT*)potential, keyword, &eps, WITH_UNITS, 1, "0.0","e", NULL);
      sprintf(keyword,  "eps_%s-%s",speciesB->name,speciesA->name);
      int flag2=object_testforkeyword((OBJECT*) potential,keyword);
      if (flag2) object_get_with_units((OBJECT*)potential, keyword, &eps, WITH_UNITS, 1, "0.0","e", NULL);
      if ((flag1+flag2)==0  && (a==b)) eps=1.0; 
      if ((flag1+flag2)==0  && (a!=b)) eps=sqrt(ljParms[aa]->eps*ljParms[bb]->eps); 

      sprintf(keyword,  "sigma_%s-%s",speciesA->name,speciesB->name);
      flag1=object_testforkeyword((OBJECT*) potential,keyword);
      if (flag1) object_get_with_units((OBJECT*)potential, keyword, &sigma, WITH_UNITS, 1, "0.0","e", NULL);
      sprintf(keyword,  "sigma_%s-%s",speciesB->name,speciesA->name);
      flag2=object_testforkeyword((OBJECT*) potential,keyword);
      if (flag2) object_get_with_units((OBJECT*)potential, keyword, &sigma, WITH_UNITS, 1, "0.0","e", NULL);
      if ((flag1+flag2)==0  && (a==b)) sigma=1.0; 
      if ((flag1+flag2)==0  && (a!=b)) sigma=0.5*(ljParms[aa]->sigma+ljParms[bb]->sigma); 

      //LJ_PARMS *ljParms_ab = ddcMalloc(sizeof(LJ_PARMS));
      LJ_PARMS *ljParms_ab = buffer+k;
      ljParms_ab->eps= eps ; 
      ljParms_ab->sigma= sigma ; 
      ljParms_ab->rcut= sigma*cutoff ; 
      LennardJones_setShift(ljParms_ab);
      ljParms[ab]  = ljParms[ba]= ljParms_ab ; 
      parms->rCutoff[ab] = parms->rCutoff[ba] = ljParms_ab->rcut; 

      if (parms->rmax < ljParms_ab->rcut)  parms->rmax = ljParms_ab->rcut; 
      if (getRank(0) == 0)
      {

         // printf("CHECK: ( %s, %s) sigma = %f eps = %f ljparms->shift = %f \n", speciesA->name, speciesB->name, sigma, eps, ljParms_ab->shift );
         //printf("eps_%s-%s = %le\n",speciesA->name,speciesB->name,ljParms_ab->eps);
         //printf("sigma_%s-%s = %le\n",speciesA->name,speciesB->name,ljParms_ab->sigma);
      }
   }
   return ljParms;
}

PAIR_PARMS *pair_parms(POTENTIAL *potential)
{
   PAIR_PARMS* parms=ddcMalloc(sizeof(PAIR_PARMS)) ;
   parms->potential=potential; 
   SYSTEM *sys= (SYSTEM*)potential->parent; 
   parms->nspecies = sys->nspecies; 
   parms->rCutoff = ddcMalloc(sizeof(double)*sys->nspecies*sys->nspecies);
   char *function; 
   object_get((OBJECT *) potential, "function", &function, STRING, 1, "NONE");
   if (strcasecmp(function,"lennardjones")==0) parms->parms = (void * )LennardJones_parms(parms); 
   if (strcasecmp(function,"lj")==0)           parms->parms = (void * )LennardJones_parms(parms); 
   if (strcasecmp(function,"TableFunction")==0) parms->parms = (void *)table_parms(potential); 
if (strcasecmp(function, "TEMPLATE_LJ_NAIVE")==0) 
{   
   printf("using tempalted strategy \n\n");
   parms->parms = (void * )LennardJones_parms(parms); 
   potential->eval_potential =(void (*) (void *,void*,void*)) pairProcessTemplatedLJ;
   potential->call_fsumX = 0;
} 
if (strcasecmp(function, "CUDA_LJ_NAIVE")==0) 
{   
   printf("using cuda lj naive strategy \n\n");
   parms->parms = (void * )LennardJones_parms(parms); 
   GPUCODE(potential->eval_potential =(void (*) (void *,void*,void*)) pairProcessTemplatedGpu;)
   potential->call_fsumX = 0;
}
if (strcasecmp(function, "CUDA_NLIST")==0) 
{   
   printf("using cuda nlist strategy \n\n");
   parms->parms = (void * )LennardJones_parms(parms); 
   GPUCODE(potential->eval_potential =(void (*) (void *,void*,void*)) pairProcessNListGpu;)
   potential->call_fsumX = 0;
   potential->use_gpu_list=1;
   potential->neighborTableType =NEIGHBORTABLE_GPU;
   potential->itype = GPU_PAIR;
}

parms->rcut= ddcMalloc(sizeof(RCUT_TYPE)*(parms->nspecies*parms->nspecies+1)); 
return parms;
}
void pair_write_dynamic(POTENTIAL *potential, FILE *file)
{
}
RCUT_TYPE *pairCutoff(SYSTEM*sys, PAIR_PARMS*parms, int *n)
{
   int ncut = 1;
   /*
      int k=0;
      for (int i = 0; i < sys->nspecies; i++) for (j = 0; j < sys->nspecies; j++) 
      {
      double rcut = parms->rcut[j+sys->nspecies*i].value ;
      parms->rcut[k].value = rcut ;
      parms->rcut[k].mode = RCUT_ALL;
      parms->rcut[k].type = i*sys->nspecies+j; 
      k++; 
      }
      */
   parms->rcut[0].value = parms->rmax;
   parms->rcut[0].mode = RCUT_ALL;
   parms->rcut[0].type = -1;
   parms->rcut[ncut].value = parms->rmax;
   parms->rcut[ncut].mode = RCUT_ALL;
   parms->rcut[ncut].type = -1;
   *n = 1; 
   return parms->rcut;
}

double getljrcuts(LJ_PARMS *parms)
{
   double rcut = parms->rcut;
   return rcut;
};

void pair(SYSTEM*sys, PAIR_PARMS *parms, ETYPE *e)
{
   int  IndexRmax = -1;
   for (int i = 0; i < sys->neighbor->nc; i++)
   {
      if (fabs(sys->neighbor->rcut[i].value - parms->rmax) < 1e-8) IndexRmax = i;
   }
   if (IndexRmax == -1 )
   {
      error_action("pair_potential error: Able to find parm rmax in neigbhor cutoff list", ERROR_IN("pair_potential", ABORT));
   }
   NPARTICLE *particles = sys->neighbor->particles;
   //TODO: This should only be called if calling a lennard jones potential
   double er=0.0;
   for (unsigned i = 0; i < sys->nlocal; i++)
   {
      int si = sys->collection->state->species[i]->index; 
      //char* namei = sys->collection->state->species[i]->name; 
      PAIRS *pij = particles[i].ifirst[IndexRmax];
      while (pij  != NULL)
      {
         RPAIR *p = pij->p;
         double r = p->r; 
         int j=pij->j;
         int sj = sys->collection->state->species[j]->index; 
         int sij = sj+sys->nspecies*si; 
         double rcut = parms->rCutoff[sij]; 
         if (r < rcut )
         {
            double dvdr; 
            double vij =  parms->fcn(parms->parms[sij],r,&dvdr) ;
            dvdr /= r; 
            er += vij;
            p->e += vij; 
            p->fp.x -= dvdr*p->x;
            p->fp.y -= dvdr*p->y;
            p->fp.z -= dvdr*p->z;
         }
         pij = pij->ilink;
      }
   }
   e->eion += er; 
}
void pairqq0(SYSTEM*sys, PAIR_PARMS *parms, ETYPE *e)
{
   double er,r,vij,dvdr,rcut; 
   NPARTICLE *particles; 
   PAIRS *pij;
   RPAIR *p; 
   int  IndexRmax = -1;
   for (int i = 0; i < sys->neighbor->nc; i++)
   {
      if (fabs(sys->neighbor->rcut[i].value - parms->rmax) < 1e-8) IndexRmax = i;
   }
   if (IndexRmax == -1 )
   {
      error_action("pair_potential error: Able to find parm rmax in neigbhor cutoff list", ERROR_IN("pair_potential", ABORT));
   }
   particles = sys->neighbor->particles;
   er=0.0;
   rcut = parms->rcut[0].value ;
   double *q = sys->collection->state->q; 
   for (unsigned i = 0; i < sys->nlocal; i++)
   {
      int si = sys->collection->state->species[i]->index; 
      pij = particles[i].ifirst[IndexRmax];
      while (pij  != NULL)
      {
         p = pij->p;
         r = p->r; 
         int j=pij->j;
         int sj = sys->collection->state->species[j]->index; 
         int sij =sj+sys->nspecies*si;
         if (r < rcut )
         {
            double qij = q[i]*q[j]; 
            vij =  qij*parms->fcn(parms->parms[sij],r,&dvdr) ;
            dvdr *= (qij/r);   
            er += vij;
            p->e += vij; 
            p->fp.x -= dvdr*p->x;
            p->fp.y -= dvdr*p->y;
            p->fp.z -= dvdr*p->z;
         }
         pij = pij->ilink;
      }
   }
   e->eion += er; 
}
/*
*/


/* Local Variables: */
/* tab-width: 3 */
/* End: */
