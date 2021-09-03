#include  <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include  <string.h>
#include <unistd.h>
#include "three_algebra.h"
#include "object.h"
#include "pio.h"
#include "ddc.h"
#include  "utilities.h"
#include "crc32.h"
#include "simulate.h"
#include "system.h"
#include "expandbuffer.h"
#include "ddcMalloc.h"
#include "io.h"
#include "box.h"
#include "units.h"
#include "state.h"
#include "mpiUtils.h"
#include "particle.h"


typedef struct velocityAutocorrelation_st {
   char *filename; 
   int last,nsample,length;
   double *vaf0,*vaf_; 
   double *msd0,*msd_; 
}  VELOCITYAUTOCORRELATION_PARMS ; 
typedef struct  ref_str {THREE_VECTOR r,v;} REF ; 

static REF *v0=NULL; 
REF *vaf_v0(void)  { return v0;}
void registerREF(REF **ref)
{
	static MPI_Datatype REF_TYPE; 
   int i,n ; 
	MPI_Aint disp[3];
	int blkcnt[] = { 6 };
	MPI_Datatype types[] = { MPI_DOUBLE };
	n = 1;
 	*ref = malloc(1*sizeof(REF));
   MPI_Get_address(&((*ref)->r), disp + 0);
   for (i = 1; i < n; i++) disp[i] -= disp[0];
   disp[0] = 0;
   MPI_Type_create_struct(n, blkcnt, disp, types, &REF_TYPE);
   MPI_Type_commit(&REF_TYPE);
   particleRegisterinfo((void *)ref, sizeof(REF), sizeof(REF), REF_TYPE,NULL,NULL);
}

VELOCITYAUTOCORRELATION_PARMS *velocityAutocorrelation_parms(ANALYSIS *analysis)
{
	VELOCITYAUTOCORRELATION_PARMS *parms ; 
	int k;
	
	parms = analysis->parms= (VELOCITYAUTOCORRELATION_PARMS *) ddcMalloc(sizeof(VELOCITYAUTOCORRELATION_PARMS));
	SIMULATE *simulate =(SIMULATE *)analysis->parent; 
	SYSTEM *sys=simulate->system; 
	object_get((OBJECT *) analysis, "filename", &parms->filename, STRING, 1, "vaf.dat");
	object_get((OBJECT *) analysis, "length", &parms->length, INT, 1, "1");
	unsigned nspecies=system_getNspecies(NULL);
	unsigned ngroups= system_getNgroup(NULL);
   if ( nspecies == 1 ) nspecies = 0;
   if ( ngroups == 1 ) ngroups = 0;
	int nv = nspecies + ngroups + 1;
	int totalLength = nv*(parms->length+1); 

// vaf0,msd0 contain the velocity autocorrelation function and mean-squared
//           displacements for the current sample.
// vaf_,msd_ contain the velocity autocorrelation function and mean-squared
//           displacements accumulated across all of the samples.
// increasing size for addition species-specific entries
	parms->vaf_ = ddcMalloc(sizeof(double)*totalLength); 
	parms->vaf0 = ddcMalloc(sizeof(double)*totalLength); 
	parms->msd_ = ddcMalloc(sizeof(double)*totalLength); 
	parms->msd0 = ddcMalloc(sizeof(double)*totalLength); 
	for (k=0;k<totalLength;k++) parms->vaf_[k] = parms->msd_[k]=0.0; 
	for (k=0;k<totalLength;k++) parms->vaf0[k] = parms->msd0[k]=0.0; 

	for (unsigned ii=0; ii<ngroups; ++ii)
	{
	   int natomsInGroup = sys->group[ii]->nMember;
	   if ( getRank(0) == 0 ) 
	   {
         printf("group = %d ; natomsInGroup = %d \n", ii, natomsInGroup);
	   }
	}

	for (unsigned ii=0; ii<nspecies; ++ii)
	{
	   int natomsInSpecies = sys->species[ii]->nMember;
	   if ( getRank(0) == 0 ) 
	   {
         printf("species = %d ; natomsInSpecies = %d \n", ii, natomsInSpecies);
	   }
	}

	parms->last=0; 
	parms->nsample=0; 
	registerREF(&v0); 
	return parms; 
}
int gcd( int a, int b )
{
	int c;
// swap a and b if a > b
	if ( a > b )
	{
      c = a; a = b; b = c;
	}
	while ( a != 0 ) 
	{
	   c = a; a = b%a;  b = c;
	}
	return b;
}
void velocityAutocorrelation_eval(ANALYSIS *analysis)
{
	
	SIMULATE *simulate; 
	SYSTEM *sys; 
	STATE *state; 
	VELOCITYAUTOCORRELATION_PARMS *parms ; 
	int i,k,nlocal; 
	parms = (VELOCITYAUTOCORRELATION_PARMS *)analysis->parms; 
	simulate =(SIMULATE *)analysis->parent; 
	sys = simulate->system; 
	state = sys->collection->state; 
	nlocal = sys->nlocal; 
	if (v0 == NULL) return ; 
// get the current velocity for each atom.
	THREE_VECTOR *__v = state_getVPtr(state) ;
// note that parms->last is really the current vaf entry, not the previous one.
	SPECIES** species = sys->collection->state->species;
	GROUP** group = sys->collection->state->group;

	unsigned nspecies=system_getNspecies(NULL);
	unsigned ngroups= system_getNgroup(NULL);
   if ( nspecies == 1 ) nspecies = 0;
   if ( ngroups == 1 ) ngroups = 0;
	int nv = nspecies + ngroups + 1;
   int totalLength = nv*(parms->length+1); 
	k=parms->last; 
	if ( k > 0 ) 
	{
      int vafoffset = 0;
      for (i=0;i<nv;i++)
      {
		   parms->msd0[k+vafoffset] = parms->vaf0[k+vafoffset] = 0.0; 
		   vafoffset += parms->length+1;
      }
		for (i=0;i<nlocal;i++) 
		{
// d is the distance from the original position of the atom
//    which is incremented in the integrator (according to Dave)
//    see defn of REF above
			THREE_VECTOR d = v0[i].r;
			double d2 = DOT(d,d);
			double v2 = DOT(v0[i].v,__v[i]);
// first add to total corr fns
			parms->msd0[k] += d2;
			parms->vaf0[k] += v2;
// now add to group corr fns
			if (ngroups != 0) 
			{
			   unsigned iGroup = group[i]->index;
			   vafoffset = (int)(1+iGroup)*(parms->length+1);
			   parms->msd0[k+vafoffset] += d2;
			   parms->vaf0[k+vafoffset] += v2;
			}
// now add to species corr fns
			if (nspecies != 0) 
			{
			   unsigned iSpecies = species[i]->index;
			   vafoffset = (int)(1+ngroups+iSpecies)*(parms->length+1);
			   parms->msd0[k+vafoffset] += d2;
			   parms->vaf0[k+vafoffset] += v2;
			}
		}
	}
	if (k == parms->length) 
	{
		parms->nsample++; 
		for (int l=0;l<totalLength;l++)
		{
			parms->msd_[l] += parms->msd0[l]; 
			parms->vaf_[l] += parms->vaf0[l]; 
			parms->msd0[l] = 0.0;
			parms->vaf0[l] = 0.0;
		}
		k=parms->last =0; 
	}
	if ( k == 0) 
	{
		for (i=0;i<nlocal;i++) 
		{
			v0[i].r=vzero  ;
			v0[i].v = __v[i]; 
		}
      int vafoffset = 0;
      for (i=0;i<nv;i++)
      {
		   parms->msd0[k+vafoffset] = parms->vaf0[k+vafoffset] = 0.0; 
		   vafoffset += parms->length+1;
		}
		for (i=0;i<nlocal;i++)
      {
         double v2 = DOT(__v[i],__v[i]);
// first total corr fns
         parms->vaf0[k] += v2;
// now group corr fns
			if (ngroups != 0) 
			{
			   unsigned iGroup = group[i]->index;
			   int vafoffset = (int)(1+iGroup)*(parms->length+1);
			   parms->vaf0[k+vafoffset] += v2;
			}
// now species corr fns
			if (nspecies != 0) 
			{
			   unsigned iSpecies = species[i]->index;
			   int vafoffset = (int)(1+ngroups+iSpecies)*(parms->length+1);
			   parms->vaf0[k+vafoffset] += v2;
			}
      }

	}
	parms->last++; 
}
void  velocityAutocorrelation_output(ANALYSIS *analysis)
{
	VELOCITYAUTOCORRELATION_PARMS *parms; 
   int outputrate = analysis->outputrate; 
   int eval_rate = analysis->eval_rate; 
	FILE *file;
	double *vaf_,*msd_; 
	parms = (VELOCITYAUTOCORRELATION_PARMS *)analysis->parms; 
	SIMULATE *simulate =(SIMULATE *)(analysis->parent); 
	SYSTEM *sys = simulate->system; 
// The following line just checks that an appropriate amount of data 
//    has been collected to do the output, so the output is sensible.
	if (parms->nsample*parms->length*eval_rate != outputrate) return ; 
	unsigned nspecies=system_getNspecies(NULL);
	unsigned ngroups= system_getNgroup(NULL);
   if ( nspecies == 1 ) nspecies = 0;
   if ( ngroups == 1 ) ngroups = 0;
	int nv = nspecies + ngroups + 1;
   int totalLength = nv*(parms->length+1); 
//	vaf_ = parms->vaf0 ; //ddcMalloc(sizeof(double)*(parms->length+1)); 
//	msd_ = parms->msd0 ; //ddcMalloc(sizeof(double)*(parms->length+1)); 
	vaf_ = ddcMalloc(sizeof(double)*totalLength); 
	msd_ = ddcMalloc(sizeof(double)*totalLength); 
 	MPI_Reduce(parms->vaf_, vaf_, totalLength, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
 	MPI_Reduce(parms->msd_, msd_, totalLength, MPI_DOUBLE, MPI_SUM, 0, COMM_LOCAL);
	double nglobal = simulate->system->nglobal; 

	double time_convert = units_convert(1.0,NULL,"t"); 
	double v2_convert = units_convert(1.0,NULL,"velocity^2"); 
	double r2_convert = units_convert(1.0,NULL,"l^2"); 
	CreateSnapshotdir(simulate, NULL);
	if (getRank(0) == 0 ) 
	{
		char filename[1024]; 
		snprintf(filename, 1023,"%s/%s", simulate->snapshotdir,parms->filename);
		file = fopen(filename, "w");

		fprintf(file, "%-33s", "#time (fs)  System vaf MSD");
		for (unsigned ii=0; ii<ngroups; ++ii)
		{
		   GROUP* g = group_by_index(NULL, ii);
		   assert(g != NULL);
		   char temp[51];
		   sprintf(temp, "  Group %s vaf MSD", g->name);
		   fprintf(file, "%-26s", temp);
		}
		for (unsigned ii=0; ii<nspecies; ++ii)
		{
		   SPECIES* s = species_by_index(NULL, ii);
		   char temp[51];
		   sprintf(temp, "  Species %s vaf MSD", s->name);
		   fprintf(file, "%-26s", temp);
		}
		fprintf(file, " (vaf in Ang^2/fs^2; msd in Ang^2)\n");
		fflush(file);

		for (int k=0;k<=parms->length;k++)  
		{
			double time = time_convert*k*simulate->dt*eval_rate;
			fprintf(file,"%f",time); 

         int vafoffset = 0;
			double vaf  = (v2_convert*vaf_[k]/parms->nsample)/nglobal; 
			double msd  = (r2_convert*msd_[k]/parms->nsample)/nglobal; 
			fprintf(file," %e %e",vaf,msd); 
         vafoffset += parms->length+1;

			for (unsigned ii=0; ii<ngroups; ++ii)
			{
			   int natomsInGroup = sys->group[ii]->nMember;
			   vaf  = (v2_convert*vaf_[k+vafoffset]/parms->nsample)/natomsInGroup; 
			   msd  = (r2_convert*msd_[k+vafoffset]/parms->nsample)/natomsInGroup; 
			   fprintf(file," %e %e",vaf,msd); 
            vafoffset += parms->length+1;
			}

			for (unsigned ii=0; ii<nspecies; ++ii)
			{
			   int natomsInSpecies = sys->species[ii]->nMember;
			   vaf  = (v2_convert*vaf_[k+vafoffset]/parms->nsample)/natomsInSpecies; 
			   msd  = (r2_convert*msd_[k+vafoffset]/parms->nsample)/natomsInSpecies; 
			   fprintf(file," %e %e",vaf,msd); 
            vafoffset += parms->length+1;
			}

			fprintf(file,"\n"); 
		}
		fflush(file);
		fclose(file);
	}
	for (int k=0;k<totalLength;k++)
	{
      parms->vaf_[k] = parms->msd_[k] = 0.0; 
	}
	parms->nsample=0; 
	ddcFree(vaf_); 
	ddcFree(msd_); 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
