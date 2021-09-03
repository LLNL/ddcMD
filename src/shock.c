#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include "group.h"
#include "system.h"
#include "state.h"
#include "species.h"
#include "object.h"
#include "ddcMalloc.h"
#include "eq.h"
#include "io.h"
#include "gid.h"
#include "pio.h"
#include "units.h"
#include "mpiUtils.h"
#include "mpiTypes.h"
#include "check_line.h"
#include "transform.h"
#include "simulate.h"
#include "pioFixedRecordHelper.h"
#include "pioVariableRecordHelper.h"
#include "shockTransform.h"
#include "particle.h"
/*
		index = atomtype[i];
		species[i] = species_by_index(NULL, index & 0xffff);
		group[i] = group_by_index(NULL, index >> 16);
*/

/*
 *  This is a group to implement a shock drive.  It sets the velocity
 *  of the particles in the group to vx = vy = 0, vz = f(t), where f(t)
 *  is given in the input deck.
 *
 *  This group is very close in spirit to both EXTFORCE and
 *  FIXED_VELOCITY.  However, neither of those implmentations does
 *  exactly what we want here.  I don't think either of those groups is
 *  in use so I could change and/or generalized them to be a shock, but
 *  writing this group took only 10 minutes.  We can always combine them
 *  in the future if we determine that is a good idea.  
*/


typedef struct shortparticle_st {unsigned domainID, type,local_index; gid_type label; double rx, ry, rz;} SHORTPARTICLE ;
typedef struct newParticle_st { gid_type label; SPECIES *species; GROUP *group; THREE_VECTOR r;} NEWPARTICLE; 
typedef struct newMaterial_st { char *filename; gid_type gidRef; THREE_VECTOR rRef; int pbc; THREE_MATRIX h0; THREE_VECTOR corner; THREE_VECTOR reducedCorner; gid_type nglobal; int nParticles; NEWPARTICLE *particle; } NEWMATERIAL; 

typedef struct shock_parms_st
{
   SIMULATE *simulate; 
   SYSTEM *system; 
   STATE  *state; 
   FILE *file;
   int nSlabs; 
   int newIndex; 
   gid_type gidRefState ;
   gid_type gidRefNew; 
   gid_type gidOffset; 
   gid_type minLabel; 
   gid_type maxLabel; 
   SHORTPARTICLE minParticle;
   SHORTPARTICLE maxParticle;
   NEWMATERIAL newMaterial; 
   double  timeLastShift; 
   double ratioRhoEst; 
   double vParticle; 
   double vShock; 
   double nglobal; 
   double volume; 
   double L; 
   double newL; 
   double zCenter; 
   double z0;
   double z1; 
   double rhoBar; 
   double rhoBarTarget; 
   double rhoBarEstimate;
   double rzM, rzP; 
   double vzM,vzP; 
   double shift; 
   
} SHOCK_PARMS;

int compareShortParticlesByGid(const void*pA, const void*pB)
{
   SHORTPARTICLE *p1 = (SHORTPARTICLE *)pA; 
   SHORTPARTICLE *p2 = (SHORTPARTICLE *)pB; 
	if (p1->label > p2->label) return 1;
	if (p1->label < p2->label) return -1;
	return 0;
}
void shockTransform_write_dynamics(TRANSFORM *t, FILE *file)
{
   double lc = units_convert(1.0, NULL, "Angstrom");
   SHOCK_PARMS *parms=(SHOCK_PARMS *)t->parms ;
   fprintf(file,"%s %s { gidRefState=%"PRIu64"; gidRefNew=%"PRIu64"; shift=%23.15e Ang;}\n",t->name,t->objclass,parms->gidRefState,parms->gidRefNew,parms->shift*lc);
}

void zeroBin(int n, double *bin)
{
   for (int i =0;i<n;i++) bin[i] = 0.0; 
}
void reduceBin(int n, double *bin)
{
   double binLocal[n]; 
   for (int i =0;i<n;i++) binLocal[i] = bin[i]; 
   MPI_Allreduce(binLocal, bin, n, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
}
double findShift(SHOCK_PARMS *parms, double delta, int nBin, double *bin)
{
   parms->rhoBar = parms->nglobal/parms->volume; 
   double nTarget = parms->rhoBarTarget*parms->volume; 
   double nMax = parms->nglobal+bin[0]; 
   double n = nMax; 

   double shift = -delta; 
   parms->rhoBarEstimate = n/parms->volume; 
   if (nTarget >  nMax) return shift; 
   for (int i=1;i<nBin;i++) 
   {
      shift -= delta; 
      n += bin[i];
      if (n < nTarget)  
      {
         double dshift =  -(nTarget-n)*delta/bin[i];
         shift +=  dshift;
         parms->rhoBarEstimate = parms->rhoBarTarget; 
         return shift;  
      }
   }
   if (getRank(0) == 0)
   {
      printf("Didn't find a solution for shift.\n"); 
      printf("This can be fixed with a better estimate for ratioRhoEst or some more coding.\n"); 
      fflush(stdout); 
   }
   assert(0);   
   return 0.0; 
}

int shockUpdateBinState(STATE *state, double offset, double delta, double weight, int nBin, double *bin)
{
   double di = 1.0/delta; 
   int np = state->nlocal; 
   int nTotal=0; 
   for (int l=0;l<np;l++)
   {
      double z = state->rz[l]+offset; 
      double s = z * di;
      int j = s; 
      assert(j>=0 ); 
      if (j< nBin) 
      {
         bin[j] += weight;
         nTotal++; 
      }
   }
   return nTotal;
}
int shockUpdateBinNew(NEWMATERIAL new, double offset, double delta, double weight, int nBin, double *bin)
{
   double di = 1.0/delta; 
   int np = new.nParticles; 
   int nTotal=0; 
   for (int l=0;l<np;l++)
   {
      double z = new.particle[l].r.z+offset; 
      double s = z * di;
      int j = s; 
      assert( j>=0 ); 
      if (j < nBin) 
      {
         bin[j] += weight;
         nTotal++; 
      }
   }
   return nTotal; 
}
void copyParticle(STATE *state,int i,int j) 
{
   if (i == j) return; 
   state->rx[i] = state->rx[j]; 
   state->ry[i] = state->ry[j]; 
   state->rz[i] = state->rz[j]; 
   state->vx[i] = state->vx[j]; 
   state->vy[i] = state->vy[j]; 
   state->vz[i] = state->vz[j]; 
   state->fx[i] = state->fx[j]; 
   state->fy[i] = state->fy[j]; 
   state->fz[i] = state->fz[j]; 
   state->q[i]  = state->q[j]; 
   state->potentialEnergy[i]  = state->potentialEnergy[j]; 
   state->virial[i]  = state->virial[j]; 
   state->sion[i]  = state->sion[j]; 
   state->label[i]  = state->label[j]; 
   state->atomtype[i]  = state->atomtype[j]; 
   state->species[i]  = state->species[j]; 
   state->group[i]  = state->group[j]; 
   // remember the particleInfo stuff; 
}
void allocateStateForNewMaterial(int n, STATE *state)
{
   state->rx = ddcMalloc(n*sizeof(double));
   state->ry = ddcMalloc(n*sizeof(double));
   state->rz = ddcMalloc(n*sizeof(double));
   state->vx = NULL; 
   state->vy = NULL; 
   state->vz = NULL; 
   state->fx = NULL; 
   state->fy = NULL; 
   state->fz = NULL; 
   state->q = NULL; 
   state->potentialEnergy = NULL; 
   state->virial = NULL; 
   state->sion = NULL; 
   state->label = ddcMalloc(n*sizeof(gid_type));
   state->atomtype = NULL;
   state->species = ddcMalloc(n*sizeof(SPECIES *));
   state->group = ddcMalloc(n*sizeof(GROUP *));
}
void freeState(STATE *state)
{
   if (state->rx != NULL)  ddcFree(state->rx);
   if (state->ry != NULL)  ddcFree(state->ry);
   if (state->rz != NULL)  ddcFree(state->rz);
   if (state->vx != NULL)  ddcFree(state->vx);
   if (state->vy != NULL)  ddcFree(state->vy);
   if (state->vz != NULL)  ddcFree(state->vz);
   if (state->fx != NULL)  ddcFree(state->fx);
   if (state->fy != NULL)  ddcFree(state->fy);
   if (state->fz != NULL)  ddcFree(state->fz);
   if (state->q  != NULL)  ddcFree(state->q );
   if (state->potentialEnergy  != NULL)  ddcFree(state->potentialEnergy );
   if (state->virial  != NULL)  ddcFree(state->virial );
   if (state->sion  != NULL)  ddcFree(state->sion );
   if (state->label  != NULL)  ddcFree(state->label);
   if (state->atomtype  != NULL)  ddcFree(state->atomtype);
   if (state->species  != NULL)  ddcFree(state->species);
   if (state->group  != NULL)  ddcFree(state->group);
}
void refTranformNewMaterial(NEWMATERIAL *new)
{
   NEWPARTICLE *particle = new->particle; 
   double hzz = new->h0.zz; 
   double zRef = new->rRef.z; 
   int l=0; 
   for (int i = 0;i<new->nParticles; i++) 
   {
      if (l < i) particle[l]= particle[i]; 
      particle[l].r.z -= zRef; 
      if ((new->pbc & 4) != 0 && (particle[i].r.z <=  0.0)) particle[l].r.z += hzz; 
      if (particle[l].r.z >  0.0) l++; 
   }
   new->nParticles = l;
}
void copyStateToNewMaterial(int n, STATE *state, NEWMATERIAL new)
{
   NEWPARTICLE *particle = new.particle; 
   for (int i = 0;i<n; i++) 
   {
      particle[i].label= state->label[i]; 
      particle[i].species= state->species[i]; 
      particle[i].group= state->group[i]; 
      particle[i].r.x= state->rx[i]; 
      particle[i].r.y= state->ry[i]; 
      particle[i].r.z= state->rz[i]; 
   }
}
//  Compare on z value and break tie with label. 
static int compar(const void *ap, const void *bp)
{
   NEWPARTICLE *a = (NEWPARTICLE *)ap; 
   NEWPARTICLE *b = (NEWPARTICLE *)bp; 
   if (a->r.z < b->r.z ) return -1; 
   if (a->r.z > b->r.z ) return  1; 
   if (a->label < b->label ) return -1; 
   if (a->label > b->label ) return  1;    
   assert(0);  //labels should be unique, so there can be no ties. 
   return 0; 
}
void checkH0(THREE_MATRIX h0)
{
   assert(h0.zz > 0); 
   assert((SQ(h0.xx) + SQ(h0.yy) + SQ(h0.xy) + SQ(h0.yx))>0.0); 
   assert((SQ(h0.xz) + SQ(h0.yz) + SQ(h0.zy)+ SQ(h0.zx))==0.0); 
}
void  printRefState(gid_type gidRef,STATE *state )
{
   FOUR_VECTOR ref4= { 0.0, 0.0,0.0,0.0}; 
   for (int j=0;j<state->nlocal;j++)      
   {
      if ( state->label[j] == gidRef) 
      {
         ref4.x = state->rx[j]; 
         ref4.y = state->ry[j]; 
         ref4.z = state->rz[j]; 
         ref4.v += 1.0; 
         break ;
      }
   }
}
THREE_VECTOR  findGidRefState(gid_type gidRef,STATE *state )
{
   FOUR_VECTOR ref4= { 0.0, 0.0,0.0,0.0}; 
   for (int j=0;j<state->nlocal;j++)      
   {
      if ( state->label[j] == gidRef) 
      {
         ref4.x = state->rx[j]; 
         ref4.y = state->ry[j]; 
         ref4.z = state->rz[j]; 
         ref4.v += 1.0; 
         break ;
      }
   }
   FOUR_VECTOR ref4Global; 
   MPI_Allreduce(&ref4, &ref4Global, 4, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   assert(ref4Global.v == 1.0) ;   //Must find gidRef and find it only once. 
   THREE_VECTOR rRef; 
   rRef.x = ref4Global.x;
   rRef.y = ref4Global.y;
   rRef.z = ref4Global.z;
   return rRef;

}
void  findGidRefNew(gid_type gidRef,NEWMATERIAL *new )
{
   FOUR_VECTOR ref4= { 0.0, 0.0,0.0,0.0}; 
   int j; 
   for (j=0;j<new->nParticles;j++)      
   {
      if ( new->particle[j].label == gidRef) 
      {
         ref4.x = new->particle[j].r.x; 
         ref4.y = new->particle[j].r.y; 
         ref4.z = new->particle[j].r.z; 
         ref4.v += 1.0; 
         break ;
      }
   }
   FOUR_VECTOR ref4Global; 
   MPI_Allreduce(&ref4, &ref4Global, 4, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   assert(ref4Global.v == 1.0) ;   //Must find gidRef and find it only once. 
   new->rRef.x = ref4Global.x;
   new->rRef.y = ref4Global.y;
   new->rRef.z = ref4Global.z;
   new->gidRef = gidRef; 
}
NEWMATERIAL readNewMaterial(char *filename)
{

   NEWMATERIAL newMaterial; 
   newMaterial.filename = strdup(filename); 
   PFILE* pfile = Popen(filename, "r", COMM_LOCAL);
   char* datatype;
   object_get(pfile->headerObject, "nrecord", &newMaterial.nglobal,U64 , 1, "0");
   object_get(pfile->headerObject, "datatype", &datatype, STRING, 1, " ");
   object_get(pfile->headerObject, "pbc", &newMaterial.pbc, INT, 1, "7");
   object_get(pfile->headerObject, "h", &newMaterial.h0, WITH_UNITS, 9, "0 0 0 0 0 0 0 0 0 ","l",NULL);
   checkH0(newMaterial.h0); 
   object_get(pfile->headerObject, "reducedcorner", &(newMaterial.reducedCorner), DOUBLE, 3, "-0.5 -0.5 -0.5");
   newMaterial.corner = matrix_vector(newMaterial.h0, newMaterial.reducedCorner);

   unsigned nRecords=0; 
   STATE state; 
   SYSTEM *sys = system_getSystem(NULL); 
   RANDOM *random = sys->random; 
   sys->random=NULL; 
   if (strcmp(datatype, "FIXRECORDASCII") == 0)  
   {
      readline_init(pfile); 
      PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) pfile->helper;
      unsigned recordLength = helper->lrec;
      nRecords = pfile->bufsize/recordLength;
      assert(pfile->bufsize%recordLength == 0);
      allocateStateForNewMaterial(nRecords,&state);
      collection_readASCII(&state, nRecords, pfile);
   }
   if (strcmp(datatype, "VARRECORDASCII") == 0)  
   {
      readline_init(pfile);
      PIO_VARIABLE_RECORD_ASCII_HELPER* helper = (PIO_VARIABLE_RECORD_ASCII_HELPER*) pfile->helper;
      nRecords = pvrah_nRecords(pfile->buf, pfile->bufsize, helper->delimiter);
      allocateStateForNewMaterial(nRecords,&state);
      collection_readASCII(&state, nRecords, pfile);
   }
   if (strcmp(datatype, "FIXRECORDBINARY") == 0) 
   {
      PIO_FIXED_RECORD_HELPER* helper = (PIO_FIXED_RECORD_HELPER*) pfile->helper;
      unsigned recordLength = helper->lrec;
      nRecords = pfile->bufsize/recordLength;
      assert(pfile->bufsize%recordLength == 0);
      allocateStateForNewMaterial(nRecords,&state);
      collection_readBINARY(&state, recordLength, nRecords, pfile);
   }
   sys->random=random; 
   newMaterial.nParticles = nRecords;
   newMaterial.particle = ddcMalloc(nRecords*sizeof(NEWPARTICLE)); 
   copyStateToNewMaterial(nRecords,&state,newMaterial); 
   Pclose(pfile);
   freeState(&state); 
   ddcFree(datatype);
   return newMaterial; 
}
MPI_Datatype shortParticle_MPIType(void)
{
   MPI_Datatype particleType;
   SHORTPARTICLE particle;
   int n = 3;
   int blkcnt[n];
   MPI_Aint disp[n];
   MPI_Datatype types[n];
   blkcnt[0] = 3;
   blkcnt[1] = 1;
   blkcnt[2] = 3;
   MPI_Get_address(&particle.domainID, &disp[0]);
   MPI_Get_address(&particle.label, &disp[1]);
   MPI_Get_address(&particle.rx, &disp[2]);
   types[0] = MPI_UNSIGNED;
   types[1] = MPI_GID_TYPE;
   types[2] = MPI_DOUBLE;
   for (int i = n-1; i >= 0; i--)
      disp[i] -= disp[0];
   MPI_Type_create_struct(n, blkcnt, disp, types, &particleType);
   MPI_Type_commit(&particleType);
   return particleType;
}
static gid_type selectGidRefNew(NEWMATERIAL *new)
{
   NEWPARTICLE *particle = new->particle; 
   gid_type gidMinLocal = gid_max; 
   double zMinLocal =  1e300;
   double zMaxLocal = -1e300; 
   for (int j =0;j<new->nParticles;j++)
   {
      if (particle[j].label < gidMinLocal ) gidMinLocal = particle[j].label; 
      if (particle[j].r.z > zMaxLocal ) zMaxLocal = particle[j].r.z; 
      if (particle[j].r.z < zMinLocal ) zMinLocal = particle[j].r.z; 
   }
   gid_type gidMin; 
   double zMin,zMax; 
   MPI_Allreduce(&gidMinLocal, &gidMin, 1, MPI_GID_TYPE, MPI_MIN, COMM_LOCAL);
   MPI_Allreduce(&zMinLocal, &zMin, 1, MPI_DOUBLE, MPI_MIN, COMM_LOCAL);
   MPI_Allreduce(&zMaxLocal, &zMax, 1, MPI_DOUBLE, MPI_MAX, COMM_LOCAL);

   gid_type gidMaxLocal  = gidMin; 
   for (int j =0;j<new->nParticles;j++)
   {
      if (particle[j].r.z  == zMax  && particle[j].label >  gidMaxLocal  ) gidMaxLocal = particle[j].label; 
   }
   gid_type gidMax;
   MPI_Allreduce(&gidMaxLocal, &gidMax, 1, MPI_GID_TYPE, MPI_MAX, COMM_LOCAL);
   return gidMax; 
}



static void minMax(SHOCK_PARMS *parms)
{
   SYSTEM *sys = parms->system; 
   STATE *state = sys->collection->state;
   struct {double value; int rank;} zMinLocal,zMaxLocal,zMinGlobal,zMaxGlobal; 
   int nlocal = sys->nlocal; 
   gid_type minLabelLocal = gid_max; 
   gid_type maxLabelLocal = 0;
   zMinLocal.value = parms->z1; 
   zMaxLocal.value = parms->z0; 
   zMinLocal.rank = getRank(0); 
   zMaxLocal.rank = getRank(0); 
   for (int j =0;j<nlocal;j++)
   {
      if (state->label[j] < minLabelLocal) minLabelLocal = state->label[j]; 
      if (state->label[j] > maxLabelLocal) maxLabelLocal = state->label[j]; 
      if (state->rz[j] > zMaxLocal.value && state->group[j]->itype == PISTON) zMaxLocal.value = state->rz[j]; 
      if (state->rz[j] < zMinLocal.value && state->group[j]->itype == PISTON) zMinLocal.value = state->rz[j]; 
   }
   MPI_Allreduce(&minLabelLocal, &parms->minLabel, 1, MPI_GID_TYPE, MPI_MIN, COMM_LOCAL);
   MPI_Allreduce(&maxLabelLocal, &parms->maxLabel, 1, MPI_GID_TYPE, MPI_MAX, COMM_LOCAL);
   MPI_Allreduce(&zMinLocal, &zMinGlobal, 1, MPI_DOUBLE_INT, MPI_MINLOC, COMM_LOCAL);
   MPI_Allreduce(&zMaxLocal, &zMaxGlobal, 1, MPI_DOUBLE_INT, MPI_MAXLOC, COMM_LOCAL);

   MPI_Datatype shortParticleType = shortParticle_MPIType();
   unsigned domainID = getRank(0); 
   if (getRank(0) == zMinGlobal.rank) 
   {
      int index = -1; 
      for (int j =0;j<nlocal;j++)
      {
         if (state->rz[j] == zMinGlobal.value) { index = j; break;}
      }
      assert(index != -1);
      parms->minParticle.domainID = domainID; 
      parms->minParticle.type = state->atomtype[index]; 
      parms->minParticle.label = state->label[index]; 
      parms->minParticle.rx = state->rx[index]; 
      parms->minParticle.ry = state->ry[index]; 
      parms->minParticle.rz = state->rz[index]; 
   }
   MPI_Bcast(&parms->minParticle, 1, shortParticleType, zMinGlobal.rank, COMM_LOCAL);

   if (getRank(0) == zMaxGlobal.rank) 
   {
      int index = -1; 
      for (int j =0;j<nlocal;j++)
      {
         if (state->rz[j] == zMaxGlobal.value) { index = j; break;}
      }
      assert(index != -1);
      parms->maxParticle.domainID = domainID; 
      parms->maxParticle.type = state->atomtype[index]; 
      parms->maxParticle.label = state->label[index]; 
      parms->maxParticle.rx = state->rx[index]; 
      parms->maxParticle.ry = state->ry[index]; 
      parms->maxParticle.rz = state->rz[index]; 
   }
   MPI_Bcast(&parms->maxParticle, 1, shortParticleType, zMaxGlobal.rank, COMM_LOCAL);
   MPI_Type_free(&shortParticleType); 
}
double zShiftStateNew(SHOCK_PARMS *parms)
{
   int cntLocal = 0; 
   double shiftLocal=0.0; 
   for (int j=0;j<parms->newMaterial.nParticles;j++)     //This works if there is only 1 particle with z  = zmax; 
   {
      double dx = parms->maxParticle.rx - parms->newMaterial.particle[j].r.x;
      double dy = parms->maxParticle.ry - parms->newMaterial.particle[j].r.y;
      double dz = parms->maxParticle.rz - parms->newMaterial.particle[j].r.z;
      if (dx  ==  0.0 && dy  ==   0.0) 
      {
         shiftLocal = dz;
         cntLocal++;
      }
   }
   int cnt=0; 
   double shift; 
   MPI_Allreduce(&cntLocal, &cnt, 1, MPI_INT, MPI_SUM, COMM_LOCAL);
   MPI_Allreduce(&shiftLocal, &shift, 1, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   if (cnt != 1)   //Hopely this doesn't happen.  If is does for non pathelogical cases then the zShift algorithim needs to be improved. 
   {
      assert(cnt== 1); 
   }
   return shift; 
}
double zNewStateOffset(SHOCK_PARMS *parms)
{
   THREE_VECTOR refNew= { 0.0,0.0,0.0}; 
   for (int j=0;j<parms->newMaterial.nParticles;j++)      
   {
      if ( parms->newMaterial.particle[j].label == parms->gidRefNew) 
      {
         refNew.x = parms->newMaterial.particle[j].r.x; 
         refNew.y = parms->newMaterial.particle[j].r.y; 
         refNew.z = parms->newMaterial.particle[j].r.z; 
         break ;
      }
   }
   STATE *state = parms->state; 
   THREE_VECTOR refState = { 0.0,0.0,0.0}; 
   for (int j=0;j<state->nlocal;j++)     
   {
      if ( parms->state->label[j] == parms->gidRefState) 
      {
         refState.x = state->rx[j]; 
         refState.y = state->ry[j]; 
         refState.z = state->rz[j]; 
         break ;
      }
   }
   THREE_VECTOR diffLocal; 
   VOP2(diffLocal,=,refState,-,refNew); 
   THREE_VECTOR diff; 
   MPI_Allreduce(&diffLocal, &diff, 3, MPI_DOUBLE, MPI_SUM, COMM_LOCAL);
   double ratio = sqrt(SQ(diff.x) + SQ(diff.y))/parms->L; 
   assert(ratio < 1e-10) ;
   return diff.z; 
}
static void shiftState(STATE *state, double shift)
{
   for (int j =0;j<state->nlocal;j++)
   {
      state->rz[j] += shift;  
   }
}
static int  markForDeletion(STATE *state, double z0, double z1)
{
   unsigned nLocal =state->nlocal; 
   for (int j =0;j<state->nlocal;j++)
   {
      if (state->rz[j] <  z0 || state->rz[j] > z1)    
      {
		  state->label[j] = (1LL<<62)-1;
         nLocal--; 
      }
   }
   return nLocal; 
}
void copyStateToShortParticle(int n,STATE *state, SHORTPARTICLE *shortParticle) 
{
   for (int i=0;i<n;i++)
   {
      shortParticle[i].domainID   = getRank(0); 
      shortParticle[i].type   = state->atomtype[i];
      shortParticle[i].local_index  = i; 
      shortParticle[i].label = state->label[i] ;
      shortParticle[i].rx = state->rx[i];      
      shortParticle[i].ry = state->ry[i];       
      shortParticle[i].rz = state->rz[i]; 
   }
}
void copyShortParticleToState(int n,STATE *state, SHORTPARTICLE *shortParticle) 
{
   for (int i=0;i<n;i++)
   {
      state->label[i]= shortParticle[i].label   ;
      state->atomtype[i]= shortParticle[i].type   ;
		state->species[i] = species_by_index(NULL, shortParticle[i].type & 0xffff);
		state->group[i] = group_by_index(NULL, shortParticle[i].type >> 16);
      state->rx[i]= shortParticle[i].rx ;     
      state->ry[i]= shortParticle[i].ry ;       
      state->rz[i]= shortParticle[i].rz ; 
   }
}
gid_type findGidOffset(unsigned nFill,gid_type maxLabel)
{

   LONG64 gidOffset=0,nFill64=nFill;
   MPI_Scan(&nFill64, &gidOffset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, COMM_LOCAL);
   assert(maxLabel+1 > nFill);
   gidOffset += maxLabel+1-nFill; 
   return gidOffset; 
}
void selectRefPair(unsigned nFill, NEWPARTICLE *particle, gid_type gidOffset, gid_type *gidRefState, gid_type *gidRefNew)
{

   double zMax=0.0;
   double zMaxLocal = -1e300;
   if(nFill > 0) zMaxLocal = particle[nFill-1].r.z;
	MPI_Allreduce(&zMaxLocal, &zMax, 1, MPI_DOUBLE, MPI_MAX, COMM_LOCAL);

   struct {gid_type gid;int rank;} local={invalidGid,-1},global;
   *gidRefState=invalidGid; 
   int rank = getRank(0);
   for (unsigned i=0;i<nFill;i++)   //Try to find a gid for which the particle z==ZMax. If more than one choose the one with the largest label/gid
   {
      if (particle[i].r.z == zMax && (particle[i].label+1 > local.gid+1 || local.rank == -1)) 
      {
         local.gid = particle[i].label;
         local.rank = rank;
         *gidRefState= gidOffset + i; 
      }
   }
	local.gid++;
   MPI_Allreduce(&local.gid, &global.gid, 1, MPI_GID_TYPE, MPI_MAX, COMM_LOCAL); 
	local.gid--;
	global.gid--;

   local.rank = -1;  //ewd need to reset this to avoid errors in high symmetry systems
   for (unsigned i=0;i<nFill;i++) 
   {
      if (particle[i].label == global.gid) 
      {
         local.rank = rank; 
         break;
      }
   }
   MPI_Allreduce(&local.rank, &global.rank, 1, MPI_INT, MPI_MAX, COMM_LOCAL);
	if(global.rank == -1) {
	  if(rank == 0) {
		 typedef unsigned long long int ull;
		 printf("invalidGid = %llu, zMax = %.5e, global.gid = %llu\n",
				  (ull) invalidGid,
				  zMax,
				  (ull) global.gid);
		 MPI_Barrier(COMM_LOCAL);
	  }
	}
   assert(global.rank != -1);  //if gidMax.rank == -1 then some how a gidRefNew was not found; 
   *gidRefNew = global.gid; 
   MPI_Bcast(gidRefState, 1, MPI_GID_TYPE, global.rank, COMM_LOCAL);
}
static int fillBox(STATE *state, NEWMATERIAL newMaterial, double z1, double offset, gid_type maxLabel, gid_type *gidRefState, gid_type *gidRefNew)
{
   unsigned nParticles = (unsigned)newMaterial.nParticles;
   NEWPARTICLE *particle = newMaterial.particle;
   unsigned nFill = 0; 
   for (; nFill< nParticles; nFill++)  
   {
      if (particle[nFill].r.z + offset >   z1)   break;
   }
   int nFillTotal=0; 
   MPI_Allreduce(&nFill, &nFillTotal, 1, MPI_UNSIGNED, MPI_SUM, COMM_LOCAL);
   int nlocal = state->nlocal; 
   if (nFillTotal > 0) 
   {
      gid_type gidOffset = findGidOffset(nFill,maxLabel); 
      selectRefPair(nFill, particle, gidOffset, gidRefState, gidRefNew);
      resize(nlocal + nFill, 2, state);
      RANDOM *random = system_getRandom(NULL); 
      for (unsigned i=0;i<nFill;i++) 
      {
         state->label[nlocal+i]   = gidOffset + i;
         state->group[nlocal+i] =   particle[i].group;
         state->species[nlocal+i] = particle[i].species;
         state->atomtype[nlocal+i] =       state->species[i]->index + (state->group[i]->index << 16);
         state->rx[nlocal+i]      = particle[i].r.x; 
         state->ry[nlocal+i]      = particle[i].r.y; 
         state->rz[nlocal+i]      = particle[i].r.z+offset; 
         state->vx[nlocal+i]      = 0.0;
         state->vy[nlocal+i]      = 0.0;
         state->vz[nlocal+i]      = 0.0;
         if (random != NULL) 
         {
            void *randomParms = random_getParms(random, nlocal+i);
            random->defaultValue(0,state->label[nlocal+i]^random->seed, randomParms);
         }
      }
   }
   nlocal += nFill; 
   return nlocal; 
}

void* shockTransform_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   char *filename; 
   char *groupName,*type; 
   SHOCK_PARMS *parms = ddcCalloc(1, sizeof(SHOCK_PARMS));
   if (!(object_testforkeyword(obj, "gidRefState")==1) || !(object_testforkeyword(obj, "gidRefNew")==1))
   {
      WAIT(0); 
      if (getRank(0)==0) {printf("\n********Check to insure gidRefState & gidRefNew are set in the %s transform object********\n\n",obj->name); fflush(stdout); }
      sleep(3); 
      WAIT(0); 
      assert(1 == 0); 
   }

   object_get((OBJECT *) obj, "rhoBarTarget", &parms->rhoBarTarget, WITH_UNITS, 1, "0.0", "1/l^3",NULL);
   object_get((OBJECT *) obj, "newMaterial", &filename, STRING, 1, "./newMaterial/atoms#");
   object_get((OBJECT *) obj, "gidRefState", &parms->gidRefState, U64, 1, invalidGidString);
   object_get((OBJECT *) obj, "gidRefNew", &parms->gidRefNew, U64, 1, invalidGidString);
   object_get((OBJECT *) obj, "nSlabs", &parms->nSlabs, INT, 1, "4");
   object_get((OBJECT *) obj, "shift", &parms->shift, WITH_UNITS, 1, "0.0", "l",NULL);
   object_get((OBJECT *) obj, "ratioRhoEst", &parms->ratioRhoEst, WITH_UNITS, 1, "0.0", "1/Ang^3",NULL);
   object_get((OBJECT *) obj, "piston", &groupName, STRING, 1, "piston");
   OBJECT *piston=object_find(groupName, "GROUP");
   object_get(piston,"type",&type,STRING,1,"none"); 
   assert(strcmp(type,"PISTON") == 0) ;
   object_get((OBJECT *) piston, "vz", &parms->vParticle, WITH_UNITS, 1, "0.0", "l/t",NULL);
   parms->simulate =(SIMULATE *)(transform->parent); 
   parms->system = parms->simulate->system;
   parms->state = parms->system->collection->state; 
   parms->newMaterial = readNewMaterial(filename); 
   parms->newL = parms->newMaterial.h0.zz; 
   parms->timeLastShift = parms->system->time; 
   registerDynamicWrite((void*)transform,(void (*)(void *,FILE *))shockTransform_write_dynamics); 
   assert(parms->rhoBarTarget > 0); 
   if (getRank(0) == 0) 
   {
      parms->file = fopen("shock.data","a"); 
      char hstring[512];
      hstring[0] = (char)0;
      sprintf(hstring+strlen(hstring),"#%-8s","loop");
      sprintf(hstring+strlen(hstring)," %16s","time (ps)    ");
      sprintf(hstring+strlen(hstring)," %12s", "gidRefState");
      sprintf(hstring+strlen(hstring)," %12s", "gidRefNew");
      sprintf(hstring+strlen(hstring)," %12s", "nAdd");
      sprintf(hstring+strlen(hstring)," %12s", "nSub");
      sprintf(hstring+strlen(hstring)," %12s", "nGlobal");
      sprintf(hstring+strlen(hstring)," %12s", "shift (Ang)");
      sprintf(hstring+strlen(hstring)," %12s", "zShock (Ang)");
      sprintf(hstring+strlen(hstring)," %12s", "Up (Ang/fs)");
      sprintf(hstring+strlen(hstring)," %12s", "Us (Ang/fs)");
      sprintf(hstring+strlen(hstring)," %12s", "rhoT(Ang^-3)");
      sprintf(hstring+strlen(hstring)," %12s", "rhoE(Ang^-3)");
      sprintf(hstring+strlen(hstring)," %12s", "rho (Ang^-3)");
      sprintf(hstring+strlen(hstring)," %12s", "rhoA(Ang^-3)");
      sprintf(hstring+strlen(hstring)," %12s", "rhoB(Ang^-3)");
      fprintf(parms->file, "%s\n", hstring);
   }


   return parms;
}

void shockTransform(TRANSFORM* transform)
{
   double lc = units_convert(1.0, NULL, "Angstrom");
   SHOCK_PARMS* parms = transform->parms;
   SYSTEM* sys = parms->system;
   STATE* state = parms->state;
   BOX_STRUCT  *box=sys->box; 
   parms->nglobal = sys->nglobal; 
   parms->volume = box->volume; 
   parms->L = box->h0.zz; 
   double cZ = box->corner.z; 
   parms->z0 = cZ; 
   parms->z1 = cZ + parms->L; 
   NEWMATERIAL new=parms->newMaterial;  
   printRefState(parms->gidRefState,state );

   findGidRefNew(parms->gidRefNew,&new);
   refTranformNewMaterial(&new);
   qsort(new.particle,new.nParticles,sizeof(NEWPARTICLE),compar); 

   THREE_VECTOR rRefState = findGidRefState(parms->gidRefState,state); 

   THREE_VECTOR diff; 
   VOP2(diff,=,rRefState,-,new.rRef); 
   double ratio = sqrt(SQ(diff.x) + SQ(diff.y))/parms->L; 
   if (ratio >= 1e-10 && getRank(0) == 0) 
   {
      printf("%"PRId64" %e\n",sys->loop,ratio);
      printf("%"PRIu64" %f %f %f\n",parms->gidRefState,rRefState.x,rRefState.y,rRefState.z);
      printf("%"PRIu64" %f %f %f\\n",new.gidRef,new.rRef.x,new.rRef.y,new.rRef.z);
   }
   assert(ratio < 1e-10) ;

   minMax(parms);

   double shiftToFill = parms->z0-parms->minParticle.rz; 
   double dSlab =  -shiftToFill;  // 
   parms->vParticle = dSlab/(parms->simulate->dt*transform->rate); 
   double vShockEst = parms->vParticle * parms->ratioRhoEst/(parms->ratioRhoEst-1.0) ;
   double shiftEst = vShockEst*(sys->time-parms->timeLastShift); 
   int  nBin  = 4*shiftEst/dSlab + 1.0; 
   if (nBin < 10) nBin = 10; 

   double bin[nBin+1]; 
   zeroBin(nBin,bin); 
   shockUpdateBinNew(new,  0, dSlab, 1.0, nBin,bin);
   int nA = shockUpdateBinState(state,-cZ, dSlab,-1.0, nBin,bin);
   bin[nBin] = nA; 
   reduceBin(nBin+1,bin); 
   double volA = (nBin-1)*dSlab*parms->volume/parms->L;  
   double rhoA = bin[nBin]/volA; 
   THREE_MATRIX hinv; 
   double rhoB = new.nglobal/matinv(&new.h0,&hinv);
   parms->vShock = parms->vParticle * rhoA/(rhoA-rhoB);
   if (getRank(0) == 0) 
   {
      double rC = units_convert(1.0, NULL, "1/Angstrom^3");
      FILE *file =stdout; 
      double n= sys->nglobal;  
      double rho = n/parms->volume; 
      fprintf(file,"%d %f %f %f %f\n",0,(0)*dSlab*lc,n,rho*rC,parms->rhoBarTarget*rC);  fflush(file); 
      for (int i =0;i<nBin;i++) 
      {
         n += bin[i]; 
         double rho = n/parms->volume; 
         fprintf(file,"%d %f %f %f %f\n",i+1,(i+1)*dSlab*lc,n,rho*rC,parms->rhoBarTarget*rC);  fflush(file); 
      }
   }
   parms->shift = findShift(parms, dSlab, nBin, bin);
   shiftState(state, parms->shift);
   double offset = rRefState.z+parms->shift; 
   unsigned nlocalUpdated; 
   nlocalUpdated=fillBox(state, parms->newMaterial,parms->z1, offset, parms->maxLabel, &parms->gidRefState, &parms->gidRefNew);
   unsigned nAddLocal = nlocalUpdated-state->nlocal;
   state->nlocal=nlocalUpdated; 
   nlocalUpdated =markForDeletion(state, parms->z0, parms->z1);
   unsigned nSubLocal = state->nlocal-nlocalUpdated;
   SHORTPARTICLE shortParticle[state->nlocal]; 
   copyStateToShortParticle(state->nlocal,state,shortParticle); 
   qsort(shortParticle,state->nlocal,sizeof(SHORTPARTICLE),compareShortParticlesByGid);
   copyShortParticleToState(state->nlocal,state,shortParticle); 
   particleSortinfo((char *)&(shortParticle->local_index),sizeof(SHORTPARTICLE),state->nlocal);
   sys->nlocal=sys->collection->size=state->nlocal=nlocalUpdated;
   gid_type nglobal = state->nlocal; 
   assert(sizeof(nglobal) == sizeof(uint64_t));
   unsigned nAdd =0; 
   unsigned nSub =0; 
   MPI_Allreduce(&nAddLocal, &nAdd, 1, MPI_UNSIGNED, MPI_SUM, COMM_LOCAL);
   MPI_Allreduce(&nSubLocal, &nSub, 1, MPI_UNSIGNED, MPI_SUM, COMM_LOCAL);
   MPI_Allreduce(&nglobal, &sys->nglobal, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, COMM_LOCAL);
   parms->nglobal = sys->nglobal; 

   double rhoBarNew = parms->nglobal/parms->volume; 
   rRefState = findGidRefState(parms->gidRefState,state); 
   parms->timeLastShift = parms->system->time; 
   if (getRank(0) == 0)
   {
      double lC = units_convert(1.0, NULL, "Angstrom");
      double vC = units_convert(1.0, NULL, "Angstrom/fs");
      double rC = units_convert(1.0, NULL, "1/Angstrom^3");
      double tC = units_convert(1.0, NULL, "fs");
      SIGNED64 loop = sys->loop; 
      double time = sys->time*tC; 
      gid_type gidRefState=parms->gidRefState; 
      gid_type gidRefNew=parms->gidRefNew; 
      gid_type nGlobal=sys->nglobal; 
      double zShock = (parms->z0 + (rhoBarNew-rhoB)/(rhoA-rhoB) * parms->L)*lC; 
      double shift=parms->shift*lC; 
      double Up = parms->vParticle*vC; 
      double Us = parms->vShock*vC; 
      double rhoT= parms->rhoBarTarget*rC; 
      double rhoE= parms->rhoBarEstimate*rC; 
      double rho = rhoBarNew*rC; 
      rhoA *= rC; 
      rhoB *= rC; 
      fprintf(parms->file,"%12.12"PRId64" %16.6f %12"PRIu64" %12"PRIu64" %12u %12u %12"PRIu64" %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
            loop,time,gidRefState,gidRefNew,nAdd,nSub,nGlobal,shift, zShock,Up,Us,rhoT,rhoE,rho,rhoA,rhoB); 
      fflush(parms->file); 
   }
}
typedef struct append_parms_st
{
   SYSTEM *system; 
   STATE  *state; 
   gid_type maxLabel; 
   NEWMATERIAL newMaterial; 
   gid_type gidRefState;
   gid_type gidRefNew;
   double gap; 
   double hzz; 

} APPEND_PARMS;
void appendTransform_write_dynamics(TRANSFORM *transform, FILE *file)
{
   APPEND_PARMS *parms=(APPEND_PARMS *)transform->parms ;
   fprintf(file,"%s %s { gidRefState=%"PRIu64"; gidRefNew=%"PRIu64";}\n",transform->name,transform->objclass,parms->gidRefState,parms->gidRefNew); fflush(file); 
}
void *appendTransform_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   APPEND_PARMS *parms = ddcCalloc(1, sizeof(APPEND_PARMS));
   char *filename; 
   object_get((OBJECT *) obj, "filename", &filename, STRING, 1, "./append/atoms#");
   object_get((OBJECT *) obj, "gap", &parms->gap, WITH_UNITS, 1, "0.0", "l",NULL);
   object_get((OBJECT *) obj, "hzz", &parms->hzz, WITH_UNITS, 1, "-1.0", "l",NULL);
   parms->system = system_getSystem(NULL);
   parms->state = parms->system->collection->state; 
   parms->newMaterial.particle =  NULL; 
   parms->newMaterial = readNewMaterial(filename); 
   if (parms->hzz < 0.0) parms->hzz = parms->newMaterial.h0.zz;
   qsort(parms->newMaterial.particle,parms->newMaterial.nParticles,sizeof(NEWPARTICLE),compar); 
   BOX_STRUCT *box = parms->system->box; 
   int nParticle=0; 
   for (; nParticle < parms->newMaterial.nParticles; nParticle++)
   {
      double rz   = parms->newMaterial.particle[nParticle].r.z -parms->newMaterial.corner.z + box->h0.zz + parms->gap; 
      if  (rz > parms->hzz)  break; 
   }
   parms->newMaterial.nParticles =nParticle; 
   parms->gidRefNew = selectGidRefNew(&parms->newMaterial);
   findGidRefNew(parms->gidRefNew,&parms->newMaterial);
   registerDynamicWrite((void*)transform,(void (*)(void *,FILE *))appendTransform_write_dynamics); 
   return parms;
}
void appendTransform(TRANSFORM* transform)
{
   APPEND_PARMS* parms = transform->parms;
   SYSTEM* sys = parms->system;
   STATE* state = parms->state;
   RANDOM *random = sys->random; 
   BOX_STRUCT *box = parms->system->box; 
   THREE_VECTOR corner = box->corner; 
   THREE_MATRIX h = parms->newMaterial.h0; 
   double zz =    box->h0.zz; 
   //double zzNew = h.zz; 
   h.zz = zz; 
   box->hfac =  matrix_matrix(h,box->hinv); 
   if (matrix_equal_tol(box->hfac, I_3x3, 1e-14)) box->hfac=I_3x3;
   //h.zz = zz+parms->gap+zzNew; 
   //box_put(box, HO, &h);
   gid_type maxLabelLocal = 0; 
   int nlocal=state->nlocal; 
   for (int kk = 0; kk < nlocal; kk++)
   {
      THREE_VECTOR rold; 
      rold.x=state->rx[kk];
      rold.y=state->ry[kk];
      rold.z=state->rz[kk];
      THREE_VECTOR r = matrix_vector(box->hfac,rold); 
      state->rx[kk]=r.x;
      state->ry[kk]=r.y;
      state->rz[kk]=r.z-corner.z+box->corner.z;
      if (state->label[kk] > maxLabelLocal) maxLabelLocal = state->label[kk]; 

   }
   MPI_Allreduce(&maxLabelLocal, &parms->maxLabel, 1, MPI_GID_TYPE, MPI_MAX, COMM_LOCAL);
   resize(nlocal+parms->newMaterial.nParticles,2,state); 
   parms->gidRefState = parms->gidRefNew + parms->maxLabel; 
   for (int kk = 0; kk < parms->newMaterial.nParticles; kk++)
   {
      state->label[nlocal]   = parms->newMaterial.particle[kk].label + parms->maxLabel;
      state->group[nlocal]   = parms->newMaterial.particle[kk].group;
      state->species[nlocal] = parms->newMaterial.particle[kk].species;
      state->rx[nlocal]      = parms->newMaterial.particle[kk].r.x; 
      state->ry[nlocal]      = parms->newMaterial.particle[kk].r.y; 
      state->rz[nlocal]      = parms->newMaterial.particle[kk].r.z -parms->newMaterial.corner.z + box->corner.z + zz + parms->gap; 
      state->vx[nlocal]      = 0.0;
      state->vy[nlocal]      = 0.0;
      state->vz[nlocal]      = 0.0;
      if (random != NULL) 
      {
       void *randomParms = random_getParms(random, nlocal);
       random->defaultValue(0,state->label[nlocal]^random->seed, randomParms);
      }
      nlocal++;
   }
   sys->nlocal=sys->collection->size=state->nlocal=nlocal ;
   gid_type nglobal = nlocal; 
   MPI_Allreduce(&nglobal, &sys->nglobal, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, COMM_LOCAL);
   double minRzLocal=  1e300; 
   double maxRzLocal= -1e300; 
   for (int kk=0;kk<nlocal;kk++) 
   {
      if ( state->rz[kk] < minRzLocal) minRzLocal = state->rz[kk]; 
      if ( state->rz[kk] > maxRzLocal) maxRzLocal = state->rz[kk]; 
   }
   double minRz,maxRz; 
   MPI_Allreduce(&minRzLocal, &minRz, 1, MPI_DOUBLE, MPI_MIN, COMM_LOCAL);
   MPI_Allreduce(&maxRzLocal, &maxRz, 1, MPI_DOUBLE, MPI_MAX, COMM_LOCAL);
   MPI_Barrier(COMM_LOCAL); 
   h.zz = maxRz-minRz; 
   box_put(box, HO, &h);
   box_put(box, HFAC, (void*)&I_3x3);
   for (int kk=0;kk<nlocal;kk++) 
   {
      state->rz[kk] -= 0.5*(maxRz+minRz); 
      if (state->rz[kk] > 0.5*h.zz) 
      {
         assert((state->rz[kk]/h.zz - 0.5)<1e-8);
         state->rz[kk] = 0.5*h.zz; 
      }
      if (state->rz[kk] < -0.5*h.zz) 
      {
         assert((-state->rz[kk]/h.zz-0.5 )<1e-8);
         state->rz[kk] = -0.5*h.zz; 
      }
   }

}
typedef struct assignGroups_parms_st
{
   SYSTEM *system; 
   STATE  *state; 
   int nGroups; 
   GROUP **groups;
   double edges[16]; 

} ASSIGNGROUPS_PARMS;
void* assignGroupsTransform_parms(TRANSFORM* transform)
{
   OBJECT* obj = (OBJECT*) transform;
   ASSIGNGROUPS_PARMS *parms = ddcCalloc(1, sizeof(ASSIGNGROUPS_PARMS));
   parms->system = system_getSystem(NULL);
   parms->state = parms->system->collection->state; 

   int nEdges = object_get((OBJECT *) obj, "edges", parms->edges, WITH_UNITS, 16, "0.0", "l",NULL);
   int nGroups = nEdges +1; 
   char *groupNames[nGroups]; 
   parms->nGroups = object_get((OBJECT *)obj, "groups", groupNames, STRING, nGroups, "");
   assert(nGroups == parms->nGroups); 
   parms->groups = ddcMalloc(nGroups*sizeof(GROUP *));
   for (int i =0;i< nGroups;i++)
   {
      parms->groups[i]= group_find(parms->system->group,groupNames[i]); 
   }


   return parms;
}
void assignGroupsTransform(TRANSFORM *transform)
{
   ASSIGNGROUPS_PARMS* parms = transform->parms;
   STATE* state = parms->state;
   BOX_STRUCT *box = parms->system->box; 
   double z0 = box->corner.z; 
   double z1 = box->h0.zz+box->corner.z; 
   int nEdges = parms->nGroups-1; 
   double knots[nEdges+1]; 
   for (int i = 0 ; i < nEdges; i++) 
   {
      if ( parms->edges[i] >= 0) knots[i] = z0 + parms->edges[i]; 
      else  knots[i] = z1 + parms->edges[i]; 
      if (i > 0) assert(knots[i] > knots[i-1]); 
   }

   for (int kk = 0; kk < state->nlocal; kk++)
   {
      int binIndex=0; 
      double rz = state->rz[kk]; 
      assert(z0 <= state->rz[kk] && state->rz[kk] <= z1  );
      for (binIndex = 0 ; binIndex <nEdges; binIndex++ )  if (rz <=knots[binIndex]  ) break; 
      GROUP *group = parms->groups[binIndex]; 
      //if (group != state->group[kk]) 
      {
         if (group->defaultValue != NULL) group->defaultValue(group, state->label[kk], kk); 
         state->group[kk]   =  group;
         state->atomtype[kk] =       state->species[kk]->index + (state->group[kk]->index << 16);
      }
   }
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
