#include "system.h" 
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "simulate.h"
#include "species.h"
#include "ddc.h"
#include "box.h"
#include "neighbor.h"
#include "collection.h"
#include "ddcMalloc.h"
#include "mpiUtils.h"
#include "energyInfo.h"
#include  "HAVEGPU.h"

void kinetic_terms(SYSTEM *, int);
void ddckinetic(int mode, SYSTEM*system);
NBR *neighbor_init(void *,char *name);

static SYSTEM *current_system = NULL;

static char* augmentRandomName( char* randomName, char** groupNames, int nGroups);

void system_write(char *string)
{
	int i,n; 
	n=current_system->nspecies ;
   for (i = 0; i < n; i++) printf("%s: =%s\n",string,current_system->species[i]->name);
}


/** Issues of initialization order:
 *  - random must be inialized before groups.
 */
unsigned system_pCalculate(SYSTEM *sys)
{
   unsigned pCalculate =0; 
   for (int j=0;j<3;j++)
   {
      int nRate; 
      int *rate;
      unsigned mask; 
      switch (j)
      {
         case 0: 
            rate=sys->pPotentialEnergyRate+1;
            nRate = sys->pPotentialEnergyRate[0];
            mask = pPotentialEnergyMask; 
            break; 
         case 1: 
            rate=sys->pVirialRate+1;
            nRate = sys->pVirialRate[0];
            mask = pVirialMask; 
            break; 
         case 2: 
            rate=sys->pSionRate+1;
            nRate = sys->pSionRate[0];
            mask = pSionMask; 
            break; 
      }
      for (int i=0;i<nRate;i++)
      {
         if ( sys->loop%rate[i] == 0) 
         {
            pCalculate ^= mask; 
            break; 
         }
      }
   }
   sys->pCalculate = pCalculate; 
   return pCalculate; 
}
SYSTEM *system_init(void *parent,char *name)
{
   char *type, *collection,*box;
   static SYSTEM *system;

   system = (SYSTEM *) object_initialize(name, "SYSTEM", sizeof(SYSTEM));
   current_system = system;
   system->parent = parent; 

   object_get((OBJECT *) system, "type", &(type), STRING, 1, NULL);
   system->loop=simulate_getLoop((SIMULATE *)parent); 
   system->time=simulate_getTime((SIMULATE *)parent); 
   system->nthreads=simulate_getNthreads((SIMULATE *)parent); 
   system->nIdleThreads = 0;
   system->moleculeClass = NULL;
   system->name = ddcCalloc(strlen(name) + 1, sizeof(char));
   system->type = ddcCalloc(strlen(type) + 1, sizeof(char));
   strcpy(system->name, name);
   strcpy(system->type, type);
   if (strcmp(type, "NORMAL") == 0)
   {
      system->itype = NORMAL;
      system->parms = ddcCalloc(1, sizeof(SYSTEM_NORMAL_PARMS));
   }

   // get names of species, groups, potentials, and randoms
   char *moleculeClass;
   char** speciesNames;
   char** groupNames;
   char** potentialNames;
   system->ngroup     = object_getv((OBJECT *) system, "groups",    (void *)&groupNames,     STRING, ABORT_IF_NOT_FOUND);
   system->npotential = object_getv((OBJECT *) system, "potential", (void *)&potentialNames, STRING, ABORT_IF_NOT_FOUND);
   int nelements;

   system->pPotentialEnergyRate = ddcMalloc(17*sizeof(int)); 
   nelements = object_get((OBJECT *) system, "pPotentialEnergyRate",    system->pPotentialEnergyRate+1,   INT, 16, "1");
   system->pPotentialEnergyRate[0] = nelements; 

   system->pVirialRate = ddcMalloc(17*sizeof(int)); 
   nelements = object_get((OBJECT *) system, "pVirialRate",    system->pVirialRate+1,   INT, 16, "1");
   system->pVirialRate[0] = nelements; 

   system->pSionRate = ddcMalloc(17*sizeof(int)); 
   nelements = object_get((OBJECT *) system, "pSionRate",    system->pSionRate+1,   INT, 16, "1");
   system->pSionRate[0] = nelements; 


   // initialization order:
   // random, species, group, box, potential, collection, neighbor   

   char* randomName;
   object_get((OBJECT *) system, "random",    &randomName,    STRING, 1, "NotSet");
   system->random=NULL; 
   randomName = augmentRandomName(randomName, groupNames, system->ngroup);
   if (strcmp(randomName,"NotSet") != 0) 
   {
      system->random = random_init(system, randomName);
      ddcFree(randomName);
   }
   object_get((OBJECT *) system, "nConstraints",   &system->nConstraints, INT, 1, "0");
   object_get((OBJECT *) system, "moleculeClass",   &moleculeClass,   STRING, 1, "NONE");
    system->energyInfo.virialCorrection = NULL;
   if (strcmp(moleculeClass,"NONE") !=0) 
   {  
      system->moleculeClass=moleculeClassInit(system,moleculeClass); 
      system->species=system->moleculeClass->species;
      system->nspecies=system->moleculeClass->nSpecies;
      system->energyInfo.virialCorrection = ddcMalloc(sizeof(THREE_VECTOR));
   }
   else
   {
      system->nspecies   = object_getv((OBJECT *) system, "species",   (void *)&speciesNames,   STRING, ABORT_IF_NOT_FOUND);
      system->species = ddcCalloc(system->nspecies + 1, sizeof(SPECIES *));
      for (int i = 0; i < system->nspecies; i++)
      {
         system->species[i] = species_init(system,speciesNames[i]);
         ddcFree(speciesNames[i]);
      }
      system->species[system->nspecies] = NULL;
      ddcFree(speciesNames);
   }
   species_put(NULL, SPECIESLIST, (void **)system->species);

   system->group = ddcCalloc(system->ngroup, sizeof(GROUP *));
   for (int i = 0; i < system->ngroup; i++)
   {
      system->group[i] = group_init(system,groupNames[i]);
      ddcFree(groupNames[i]);
   }
   ddcFree(groupNames);

   object_get((OBJECT *) system, "collection", &(collection), STRING, 1, NULL);
   object_get((OBJECT *) system, "box", &(box), STRING, 1, NULL);
   system->changed = 7;
   system->box = box_init(system,box);

   system->collection = collection_init(system,collection);
   system->nlocal = system->collection->size;
   system->nremote = 0;
   system->nion = system->nlocal + system->nremote;
   box_put(system->box,BOX_TIMESTART,&system->time);
   gid_type nLocal = system->nlocal;
   MPI_Allreduce(&nLocal, &system->nglobal, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);

   GPUCODE(allocSendGPUState(system->collection, system->collection->size);)
   GPUCODE(allocGPUBoxInfo(system);)

      if (system->nglobal != system->collection->gsize )
      {
         if (getRank(0) == 0)
         {
            printf("global=%"PRIu64" %"PRIu64"\n",system->nglobal,system->collection->gsize); 
            printf("FATAL ERROR: Number of particles read does not match number in file.\n");
         }
         WAIT(0); 
         abortAll(71);
      }
   char *energyInfoObj; 
   object_get((OBJECT *) system, "energyInfo", &(energyInfoObj), STRING, 1, "energyInfo");
   system->energyInfoObj = energyInfo_init(system,energyInfoObj);
   zeroEType(&system->energyInfo);

   system->potential = ddcCalloc(system->npotential, sizeof(POTENTIAL *));
   system->neighborTableType = 0; 
   for (int i = 0; i < system->npotential; i++)
   {
      system->potential[i] = potential_init(system, potentialNames[i]);
      system->neighborTableType |= system->potential[i]->neighborTableType; 

      ddcFree(potentialNames[i]);
   }
   ddcFree(potentialNames);

   object_get((OBJECT *) system, "neighbor", &(name), STRING, 1, NULL);
   system->neighbor = neighbor_init(system,name);
   object_get((OBJECT *) system, "volume_per_atom", &system->energyInfo.vol_per_atom, DOUBLE, 1, "-1");
   system->kineticFunc=(void(*)(void*, int))kinetic_terms;
   return system;
}

RANDOM *system_getRandom(SYSTEM*system) 
{	
   if (system==NULL) system=current_system; 
   return system->random; 
}
SIGNED64  system_getLoop(SYSTEM*system) 
{	
   if (system==NULL) system=current_system; 
   return system->loop; 
}
double system_getTime(SYSTEM*system) 
{	
   if (system==NULL) system=current_system; 
   return system->time; 
}
gid_type system_getNglobal(SYSTEM *system) 
{ 
   if (system==NULL) system=current_system; 
   return system->nglobal;
} 
unsigned system_getNlocal(SYSTEM *system) 
{ 
   if (system==NULL) system=current_system; 
   return system->nlocal;
} 
STATE *system_getState(SYSTEM *system) 
{
   if (system==NULL) system=current_system; 
   return system->collection->state; 
}
SYSTEM *system_getSystem(SYSTEM *system) 
{
   if (system==NULL) system=current_system; 
   return system; 
}
int system_getNspecies(SYSTEM *system) 
{
   if (system==NULL) system=current_system; 
   return system->nspecies; 
}
SPECIES **system_getSpecies(SYSTEM*system) 
{	
   if (system==NULL) system=current_system; 
   return system->species; 
}
int system_getNgroup(SYSTEM *system) 
{
   if (system==NULL) system=current_system; 
   return system->ngroup; 
}
GROUP **system_getGroup(SYSTEM *system) 
{
   if (system==NULL) system=current_system; 
   return system->group; 
}
int system_get(SYSTEM*system, int get, void *ptr)
{
   if (system==NULL) system=current_system; 
   switch (get)
   {
      case SYSTEMNAME:
         *(char **)ptr = system->name;
         return 1;
      case SYSTEMITYPE:
         *((enum SYSTEM_CLASS *)ptr) = system->itype;
         return 1;
      case SYSTEMTYPE:
         *(char **)ptr = system->type;
         return 1;
      case RETURN_CODE:
         *((int *)ptr) = system->return_code;
         return 1;
      case ENERGY:
         *((double *)ptr) = system->energy;
         return 1;
      case PRESSURE:
         *((double *)ptr) = system->p;
         return 1;
      case NLOCAL:
         *((unsigned *)ptr) = system->nlocal;
         return 1;
      case NGLOBAL:
         *((gid_type *)ptr) = system->nglobal;
         return 1;
      case CURRENT_SYSTEM:
         *(SYSTEM **) ptr = current_system;
         return 1;
      default:
         return 0;
   }
}

int system_put(SYSTEM*system, int put, void *ptr)
{
   if (system==NULL) system=current_system; 
   switch (put)
   {
      case SYSTEM_LOOP:
         system->loop = *(SIGNED64*)ptr; 
         return 1; 
      case SYSTEM_TIME:
         system->time = *(double*)ptr; 
         box_put(system->box,BOX_TIME,(void *)&system->time); 
         return 1; 
      default:
         return 0;
   }
}
void system_putNlocal(SYSTEM *system,unsigned nlocal) 
{ 
   if (system==NULL) system=current_system; 
   system->nlocal=nlocal;
} 
void print_atom_by_label(LONG64 atom_label,char *string)
{
   double *rx,*ry,*rz,*vx,*vy,*vz; 	
   LONG64 *label; 
   STATE *state; 
   SYSTEM *sys; 
   state = system_getState(NULL);  
   sys = system_getSystem(NULL);  
   rx = state->rx; ry = state->ry; rz = state->rz;
   vx = state->vx; vy = state->vy; vz = state->vz;
   label = state->label; 
   for (unsigned kk=0;kk<sys->nlocal;kk++) 
   {
      if (atom_label == label[kk]) 
      {
         fprintf(ddcfile,"%16s %6d: %8"PRIu64" %16"PRIu64" %8d %8d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",string,getRank(0),sys->loop,label[kk],kk,sys->nlocal,rx[kk],ry[kk],rz[kk],vx[kk],vy[kk],vz[kk]); 
         fflush(ddcfile); 
         return; 
      }
   }
}

/*
   RANDOM* system_getRandomByName(SYSTEM* sys, const char* name)
   {
   for (int ii=0; ii<sys->nrandom; ++ii)
   if (strcmp(name, sys->random[ii]->name) == 0)
   return sys->random[ii];
   return NULL;
   }
   */


/** Searches for name in list.  Returns index where found or nList if
 * not found. */
/*
   static unsigned findName(char* name, char** list, unsigned nList)
   {
   for (unsigned ii=0; ii<nList; ++ii)
   if (strcmp(name, list[ii]) == 0)
   return ii;
   return nList;
   }
   */

/** This function is strictly for backwards compatibility.
 *
 *  In April 2010 we realized that the only sensible way for multiple
 *  groups to use the same random object was for the random object to be
 *  owned at a higher level than the groups.  System is the reasonable
 *  place since it already owns similar structures like the species and
 *  group lists.  We decided that users should list the names of the
 *  RANDOM objects needed by the simulation in the random keyword of the
 *  SYSTEM object.  The system_init function would then initialize the
 *  necessary RANDOM objects.  The one problem with this scheme is that
 *  it breaks every existing ddcMD input deck.  To provide backwards
 *  compatibility we decided that before initializing RANDOM objects we
 *  will scan through all of the GROUP names given to system.  We will
 *  find any names of random objects mentioned by the groups and add
 *  those names to the list of random objects.
 *
 *  We do this for compatibility only.  We will not extend this kind of
 *  treatment to other types of objects because along that road lies
 *  madness.
 *
 *  This function scans the named group objects to find any mentioned
 *  random names and adds them to the randomNames.  We also ensure that
 *  names appear only once in the list.
 */
static int groupsWithRandom[]= {LANGEVIN, RELATIVISTICLANGEVIN, RADIATION, IONIZATION, TNBURN,-1}; 

char* augmentRandomName(char* randomName, char** groupNames, int nGroups)
{
   for (int ii=0; ii<nGroups; ++ii)
   {
      OBJECT *group = object_find(groupNames[ii], "GROUP");
      int itype = ((GROUP *)group)->itype; 
      int *grouptype = groupsWithRandom; 
      while( *grouptype++ != -1) 
      {  
         if (itype == *grouptype ) break;
         grouptype++; 
      }
      if (*grouptype != -1) 
      {
         if (object_testforkeyword(group, "random") == 1) 
         {
            char* name;
            object_get(group, "random", &name, STRING, 1, " ");
            if (strcmp(randomName,"NotSet")==0) 
            {
               free(randomName);
               randomName = strdup(name); 
            }
            if (strcmp(name,randomName) != 0) 
            {
               if (getRank(0) == 0) 
               {
                  printf("ERROR\n"); 
                  printf("The currrent code only allows for one per particle random number generator. This should be defined in the SYSTEM object.\n");
                  printf("For backwards compatibilty the random object will be scanned for in the LANGEVIN, RELATIVISTICLANGEVIN, RADIATION, IONIZATION and TNBURN groups.\n"); 
                  printf("There is a constraint that each random keyword if define must have the same value as defined in any other  of these groups.\n") ;
               }
            }
            ddcFree(name);
         }
      }
   }
   return randomName;
}

void system_getDomainParticlesEnds(SYSTEM *system, double* rmin, double* rmax)
{
   double* rx = system->collection->state->rx; 
   double* ry = system->collection->state->ry; 
   double* rz = system->collection->state->rz; 
   int nlocal = system->nlocal; 

   for(int i=0;i<3;i++)rmin[i] = 1000000.;
   for(int i=0;i<3;i++)rmax[i] =-1000000.;

   for (int ip = 0; ip < nlocal; ip++)
   {
      double x = rx[ip]; 
      double y = ry[ip]; 
      double z = rz[ip]; 

      if( x>rmax[0] )rmax[0]=x;
      if( y>rmax[1] )rmax[1]=y;
      if( z>rmax[2] )rmax[2]=z;

      if( x<rmin[0] )rmin[0]=x;
      if( y<rmin[1] )rmin[1]=y;
      if( z<rmin[2] )rmin[2]=z;
   }
   return;
}

void system_getParticlesEnds(SYSTEM *system, THREE_VECTOR* rminglobal, THREE_VECTOR* rmaxglobal)
{
   double rmin[3];
   double rmax[3];
   system_getDomainParticlesEnds(system, &rmin[0],&rmax[0]);
   MPI_Allreduce(&rmin[0], &rminglobal->x, 3, MPI_DOUBLE, MPI_MIN, COMM_LOCAL);
   MPI_Allreduce(&rmax[0], &rmaxglobal->x, 3, MPI_DOUBLE, MPI_MAX, COMM_LOCAL);
}

/** We expect the caller to pass in an actual SYSTEM pointer.  No use of
 *  current_system in this function. */
enum POT_COMM_MODE system_getCommMode(SYSTEM* sys)
{
   assert(sys);
   assert(sys->npotential>0 && sys->potential[0]);
   enum POT_COMM_MODE commMode = sys->potential[0]->commMode;
   for (int ii=0; ii<sys->npotential; ++ii)
      if (sys->potential[ii]->commMode != commMode && getRank(0) == 0)
      {
         printf("FATAL ERROR: You have specified two (or more) potentials\n"
               "with different DDC comm_modes.  You'll need to contact\n"
               "a developer to help you sort this one out.\n");
         abortAll(-1);					 
      }

   return commMode;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
