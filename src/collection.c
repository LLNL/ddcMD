#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define  _LARGE_FILE
#include "collection.h"
#include "three_algebra.h"
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include "external.h"
#include "object.h"
#include "error.h"
#include "utilities.h"
#include "ddcMalloc.h"
#include "box.h"
#include "system.h"
#include "gidShuffle.h"
#include "units.h"
#include "codata.h"
#include "mpiUtils.h"
#include "particle.h"
#include "pio.h"
#include "gid.h"

#include "simulate.h"
#include "state_redist.h"



FILE* ddcfile;



// MPI_Datatype mpiatomtype;
// extern int mpiatomtype;

static void collection_buildLattice(COLLECTION* c);
static void misspelledKeywordError(void);

COLLECTION *collection_init(void *parent, char *name)
{
	static COLLECTION *collection;

	timestampBarrier("Starting Collection_Init", COMM_LOCAL);
	collection = (COLLECTION *) object_initialize(name, "COLLECTION", sizeof(COLLECTION));
	collection->parent = parent; 
	collection->size = 0;
   collection->dynamicInfoFromFileHeader=0;
	collection->state = ddcMalloc(sizeof(STATE));
	STATE* state = collection->state;
	SYSTEM *sys=system_getSystem(NULL); 
	state->rx = state->ry = state->rz = state->vx = state->vy = state->vz = state->fx = state->fy = state->fz = NULL;
	state->q = NULL;
	state->label = NULL;
	state->group = NULL;
	state->species = NULL;
	state->atomtype = NULL;
	state->potentialEnergy = NULL;
	state->virial = NULL;
	state->sion = NULL;

	state->cost = ddcMalloc(1);
	particleRegisterinfo((void *)&state->cost, sizeof(*(state->cost)), sizeof(*(state->cost)), MPI_FLOAT,NULL,NULL);

	particleRegisterinfo((void *)&state->vx, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
	particleRegisterinfo((void *)&state->vy, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
	particleRegisterinfo((void *)&state->vz, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
/*
	particleRegisterinfo((void *)&state->q, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
	particleRegisterAccum((void *)&state->fx, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
	particleRegisterAccum((void *)&state->fy, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
	particleRegisterAccum((void *)&state->fz, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
	particleRegisterAccum((void *)&state->potentialEnergy, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
	particleRegisterAccum((void *)&state->virial, sizeof(double), sizeof(double), MPI_DOUBLE,NULL,NULL);
*/
	char* type;
	object_get((OBJECT *) collection, "type", &type, STRING, 1, "READ");

   if (strncmp(type, "LATTICE", strlen("LATTICE")) == 0  )	collection_buildLattice(collection);
   if (strncmp(type, "READ", strlen("READ")) == 0)  
   {
      char *filename;
      char *dynamicInfo; 
      object_get((OBJECT *) collection, "files", &filename, LITERAL, 1, " ");
      object_get( (OBJECT*)collection, "size", &collection->gsize, U64, 1, "0"); 
      object_get( (OBJECT*)collection, "dynamicInfo", &dynamicInfo,STRING, 1, "notFromFileHeader"); 
      if (strcmp(dynamicInfo,"fromFileHeader")==0) collection->dynamicInfoFromFileHeader=1; 
      PFILE* pfile = Popen(filename, "r", COMM_LOCAL);
      collection_read(collection,pfile); 
      Pclose(pfile);
      ddcFree(filename);
   }
   ddcFree(type);

   for (int i=0;i<collection->size;i++) state->q[i] = ((ATOMTYPE_PARMS*)(state->species[i])->parm)->charge;
   RANDOM *random = system_getRandom(sys); 
   if (random != NULL)
   {
      unsigned useDefault; 
      MPI_Allreduce(&(random->useDefault), &useDefault, 1, MPI_UNSIGNED, MPI_MAX, COMM_LOCAL);
      random->useDefault=useDefault; 
      if (useDefault)
      {
         for (int i=0;i<collection->size;i++) 
         {
            void *randomParms = random_getParms(random, i);
            random->defaultValue(0, state->label[i],randomParms); 
         }
      }
   }

   /* Optionally redistribute particles for load balancing, and generate processor map / domain centers */ 
   {
      char *parttype;
      void *simulate_p;

      state->domcen_active = 0;
      state->domcen_np = 0;
      state->domcen_np = 0;
      state->domcen_vec = NULL;

      (void) simulate_get(NULL,SIMULATE_LAST,&simulate_p);
      //         printf("{%s:%d} in %s(): simulate_p = %llu\n", __FILE__,__LINE__,__func__, (unsigned long long int) simulate_p);
      object_get((OBJECT *) simulate_p, "partitioning", &parttype, STRING, 1, "NONE");
      //       printf("{%s:%d} in %s(): partitioning = %s\n", __FILE__,__LINE__,__func__, parttype);
      if(strcmp(parttype,"RECURSIVE_BISECTION") == 0 || strcmp(parttype,"recursive_bisection") == 0) 
      {
         timestampBarrier("Performing recursive bisection particle redistribution.",COMM_LOCAL);
         collection->size = state_redist(state,COMM_LOCAL,collection->size);
      }
      ddcFree(parttype);
   }

   /*
      for (int i=0;i<sys->ngroup;i++)
      {
      unsigned useDefault; 
      GROUP *group=sys->group[i]; 
      MPI_Allreduce(&(group->useDefault), &useDefault, 1, MPI_UNSIGNED, MPI_MAX, COMM_LOCAL);
      group->useDefault=useDefault; 
      }
      for (int i=0;i<collection->size;i++) 
      {
      GROUP *group=state->group[i]; 
      if (group->useDefault && group->defaultValue!=NULL)
      group->defaultValue(group, state->label[i], i); 
      }
      */
   for (int i=0;i<collection->size;i++) 
   {
      state->group[i]->nMember++;
      state->species[i]->nMember++;
   }

   int n=sys->ngroup+sys->nspecies; 
   gid_type local[n],global[n]; 
   int k=0; 
   for (int i=0;i<sys->ngroup;i++) local[k++]= sys->group[i]->nMember ; 
   for (int i=0;i<sys->nspecies;i++)local[k++]= sys->species[i]->nMember; 
   MPI_Allreduce(local, global, n, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
   k=0;  
   for (int i=0;i<sys->ngroup;i++) sys->group[i]->nMember  = global[k++]; 
   for (int i=0;i<sys->nspecies;i++)sys->species[i]->nMember= global[k++]; 


   timestampBarrier("Finished Collection_Init", COMM_LOCAL);
   return collection;
}
void collection_read(COLLECTION *collection, PFILE *pfile)
{
   if (!object_testforkeyword(pfile->headerObject, "groups")) object_replacekeyword(pfile->headerObject,"groups","group");
   if (!object_testforkeyword(pfile->headerObject, "types")) object_replacekeyword(pfile->headerObject,"types","ATOM");

   char* datatype;
   object_get(pfile->headerObject, "datatype", &datatype, STRING, 1, " ");

   if (collection->dynamicInfoFromFileHeader)
   {
      double time; 
      SIGNED64 loop; 
      object_get(pfile->headerObject, "time", &time, WITH_UNITS, 1, "0.0","t",NULL);
      object_get(pfile->headerObject, "loop", &loop, U64, 1, "0");

      simulate_put(NULL,SIMULATE_LOOP,&loop);
      simulate_put(NULL,SIMULATE_TIME,&time);

      THREE_MATRIX h; 
      object_get(pfile->headerObject, "h", &h, WITH_UNITS, 9, "1 0 0 0 1 0 0 0 1","l",NULL);
      box_put(NULL,HO,&h);
   }
   //off_t p1 = ftello(pfile->readFile[0]); 
   if (strcmp(datatype, "BXYZ") == 0)            collection_readFIXRECORDBINARY(collection, pfile);
   if (strcmp(datatype, "FIXRECORDASCII") == 0)  collection_readFIXRECORDASCII(collection, pfile);
   if (strcmp(datatype, "VARRECORDASCII") == 0)  collection_readVARRECORDASCII(collection, pfile);
   if (strcmp(datatype, "FIXRECORDBINARY") == 0) collection_readFIXRECORDBINARY(collection, pfile);
   //off_t p2 = ftello(pfile->readFile[0]); 

   THREE_VECTOR reducedCorner0 = box_get_reducedcorner(NULL); 
   THREE_VECTOR reducedCorner; 
   char defaultString[256]; 
   sprintf(defaultString, "%18.15f %18.15f %18.15f",reducedCorner0.x,reducedCorner0.y, reducedCorner0.z);
   object_get( (OBJECT*)pfile->headerObject, "reducedCorner", &reducedCorner,DOUBLE, 3,defaultString); 
   sprintf(defaultString, "%18.15f %18.15f %18.15f",reducedCorner.x,reducedCorner.y, reducedCorner.z);
   object_get( (OBJECT*)collection, "reducedCorner", &reducedCorner,DOUBLE, 3,defaultString);   

   THREE_MATRIX h = box_get_h(NULL);
   THREE_VECTOR deltaReduced = {reducedCorner0.x - reducedCorner.x, reducedCorner0.y - reducedCorner.y,reducedCorner0.z - reducedCorner.z}; 
   THREE_VECTOR deltaCorner = matrix_vector(h, deltaReduced);
   for  (int i=0;i<collection->size;i++)
   {
      collection->state->rx[i] += deltaCorner.x; 
      collection->state->ry[i] += deltaCorner.y; 
      collection->state->rz[i] += deltaCorner.z; 
   }
   ddcFree(datatype);
}


/**
 *  This routine is constructed so that it will generate the same atomic
 *  configuration regardless of how many tasks it runs on.  This enables
 *  us to use generated configurations when testing without having to
 *  worry about getting different results for different task counts.
 *
 *  To get the same result regardless of task count we re-seed the
 *  random number generator for every atom using the lowest 48 bits of
 *  the gid XOR'd with the user provided 48 bit seed.
 */
void collection_buildLattice(COLLECTION* c)
{
   STATE* state = c->state;

   THREE_VECTOR lattice_vectors[3];
   double scale,delta,temperature; 
   char **species_list = NULL; 
   char **group_list = NULL; 
   int lattice_size[3];
   gid_type seed;
   int randomize=0;
   double vcm[3];

   int nelements;
   double aveParticlesPerRank; 
   nelements   = object_get((OBJECT *) c, "scale", &scale, WITH_UNITS, 1, "1","l",NULL);
   nelements   = object_get((OBJECT *) c, "delta", &delta, DOUBLE, 1, "0");
   nelements   = object_get((OBJECT *) c, "lattice_vectors", lattice_vectors, DOUBLE, 9, "1 0 0 0 1 0 0 0 1");
   nelements   = object_get((OBJECT *) c, "lattice_size", lattice_size, INT, 3, "1 1 1");
   int nbasis  = object_getv((OBJECT *) c, "species", (void *)&species_list, STRING,IGNORE_IF_NOT_FOUND);
   nelements   = object_get((OBJECT *) c, "temperature", &temperature, WITH_UNITS, 1, "0","K",NULL);
   nelements   = object_get((OBJECT *) c, "aveParticlesPerRank", &aveParticlesPerRank, DOUBLE, 1, "0.0");
   assert(nbasis > 0) ; 
   if (aveParticlesPerRank > 0.0) 
   {
      double size = aveParticlesPerRank*getSize(0);
      int lx = lrint(cbrt(size/nbasis)); 
      int ly = lx;
      int lz = lx;
      lattice_size[0] = lx; 
      lattice_size[1] = ly; 
      lattice_size[2] = lz; 
   }
   c->gsize = ((LONG64)nbasis) * ((LONG64)lattice_size[0]) * ((LONG64)lattice_size[1]) * ((LONG64)lattice_size[2]); 


   if (object_testforkeyword((OBJECT *) c, "group")) misspelledKeywordError();

   int nGroups = object_getv((OBJECT *) c, "groups",   (void *)&group_list, STRING,IGNORE_IF_NOT_FOUND);
   if (nGroups > 0 && nGroups != nbasis) 
   {
      error_action("Number of groups specifed is not consistent with species list in  object <<",
            c->name, ">> ", ERROR_IN("collection_buildLattice", ABORT));
   }
   THREE_VECTOR* basis = (THREE_VECTOR *)ddcMalloc((nbasis+1)*sizeof(THREE_VECTOR));
   nelements = object_get((OBJECT *) c, "seed", &seed, U64, 1, "0X1717357B31472A59");
   nelements = object_get((OBJECT *) c, "randomizeSeed", &randomize, INT, 1, "0");
   if (randomize != 0)
      seed = generateRandomSeed();
   nelements = object_get((OBJECT *) c, "basis", basis, DOUBLE, (nbasis+1)*3, "0 0 0");
   if (nelements != nbasis*3) 
   {
      error_action("Number basis elements specifed is not consistent with species list in  object <<",
            c->name, ">> ", ERROR_IN("collection_buildLattice", ABORT));
   }

   nelements = object_get((OBJECT* )c, "vcm", vcm, DOUBLE, 3, "0 0 0");
   if (nelements != 3)
      error_action("Improper vcm specified in object <<", c->name, ">>.\n"
            "Please specify three velocity components.",
            ERROR_IN("collection_buildLattice", ABORT));
   double length_convert = units_convert(1.0,"l",NULL); 
   double time_convert = units_convert(1.0,"t",NULL); 
   double velocity_convert = length_convert/time_convert;
   for (unsigned ii=0; ii<3; ++ii)
      vcm[ii] *= velocity_convert;

   if (nGroups == 0)
   {
      nGroups = nbasis;
      group_list = ddcMalloc(nbasis * sizeof(char*));
      for (int ii=0; ii<nGroups; ++ii) group_list[ii] = strdup("FREE");
   }

   THREE_VECTOR shift;
   shift.x = shift.y = shift.z = 0.0; 
   for (int ii=0; ii<nbasis; ++ii)
   {
      basis[ii].x *= scale ; 
      basis[ii].y *= scale ; 
      basis[ii].z *= scale ; 
      shift.x -= basis[ii].x/nbasis; 
      shift.y -= basis[ii].y/nbasis; 
      shift.z -= basis[ii].z/nbasis; 
   }
   THREE_VECTOR h0[3]; 
   for (unsigned ii=0; ii<3; ++ii)
   {
      lattice_vectors[ii].x *= scale; 
      lattice_vectors[ii].y *= scale; 
      lattice_vectors[ii].z *= scale; 
      h0[ii].x = lattice_vectors[ii].x * lattice_size[ii]; 
      h0[ii].y = lattice_vectors[ii].y * lattice_size[ii]; 
      h0[ii].z = lattice_vectors[ii].z * lattice_size[ii]; 
      shift.x += 0.5*lattice_vectors[ii].x ;
      shift.y += 0.5*lattice_vectors[ii].y ;
      shift.z += 0.5*lattice_vectors[ii].z ;
   }
   box_put(NULL,HO,&h0);
   box_put(NULL,HFAC,(void*)&I_3x3);
   THREE_VECTOR* corner;
   box_get(NULL,CORNER_PTR,&corner);

   long long int nb = nbasis;
   long long int lx = lattice_size[0];
   long long int ly = lattice_size[1];
   long long int lz = lattice_size[2];

   long long int  sizeGlobal = nb*lx*ly*lz;
   long long int  nTasks = (int)getSize(0);
   long long int myRank = (int)getRank(0);

   long long int  sizeLocal = (sizeGlobal/nTasks) + 1LL;
   long long int  iBase = sizeLocal*myRank;
   if (iBase + sizeLocal > sizeGlobal) sizeLocal = MAX(0LL, sizeGlobal - iBase);
   {
      gid_type send = sizeLocal;
      gid_type recv;
      MPI_Allreduce(&send, &recv, 1, MPI_GID_TYPE, MPI_SUM, COMM_LOCAL);
      assert(recv == (gid_type)sizeGlobal);
   }

   c->size = sizeLocal;
   resize(c->size, 2, state);

   int nx = lattice_size[0];
   int ny = lattice_size[1];

   assert( state->species!=NULL );
   assert( state->atomtype!=NULL );
   for (unsigned ii=0; ii<sizeLocal; ++ii)
   {
      gid_type gid = iBase + ii;
      state->label[ii] = gid; 
      PRAND48_STATE handle = prand48_init(gid, seed, 0x9753bcde8642llu);
      int ib = gid % nbasis;
      gid /= nbasis;
      int ix = gid % nx;
      gid /= nx;
      int iy = gid % ny;
      gid /= ny;
      int iz = gid;

      state->species[ii] = species_find(NULL, species_list[ib]);
      if( state->species[ii] == NULL )
         error_action("Cannot find species <<", species_list[ib], ">>.\n" "Please define species.", ERROR_IN("collection_buildLattice", ABORT));

      state->group[ii] = group_find(NULL, group_list[ib]);
      assert(state->group[ii] != NULL);

      state->atomtype[ii] = state->species[ii]->index + (state->group[ii]->index << 16);

      assert(state->species[ii] != NULL);
      state->rx[ii] = corner->x + shift.x + ix*lattice_vectors[0].x  + iy*lattice_vectors[1].x + iz*lattice_vectors[2].x + basis[ib].x+delta*scale*prand48(&handle); 
      state->ry[ii] = corner->y + shift.y + ix*lattice_vectors[0].y  + iy*lattice_vectors[1].y + iz*lattice_vectors[2].y + basis[ib].y+delta*scale*prand48(&handle); 
      state->rz[ii] = corner->z + shift.z + ix*lattice_vectors[0].z  + iy*lattice_vectors[1].z + iz*lattice_vectors[2].z + basis[ib].z+delta*scale*prand48(&handle); 
      THREE_VECTOR vv={0.0,0.0,0.0}; 
      if (temperature > 0.0) 
      {
         double mass = ((ATOMTYPE_PARMS *) (state->species[ii]->parm))->mass;
         double sigma = sqrt(kB*temperature/mass);
         vv = gasdev3d0(sigma, &handle); 
      }
      state->vx[ii] = vv.x+ vcm[0];
      state->vy[ii] = vv.y+ vcm[1];
      state->vz[ii] = vv.z+ vcm[2];
   }

   for (int ii=0; ii<nbasis; ++ii) ddcFree(species_list[ii]);
   ddcFree(species_list);
   for (int ii=0; ii<nGroups; ++ii) ddcFree(group_list[ii]);
   ddcFree(group_list);

   c->state->nlocal = sizeLocal;
   gidShuffle(c, seed);

   // initialize groups after the gid shuffle since label may used in generating a seed for the per-particle random state.
   RANDOM *random = system_getRandom(NULL); 
   if (random != NULL) random->useDefault = 1; 
   for (unsigned ii=0; ii<sizeLocal; ++ii) if (state->group[ii]->parse != NULL) state->group[ii]->parse(state->group[ii], NULL, ii);   //JNG CHECK 

}

void misspelledKeywordError(void)
{
   if (getRank(0) != 0)
      return;

   printf("\nIt appears that you are trying to use the 'group' keyword in the \n"
         "COLLECTION object. This was the correct spelling prior to version 1638.  \n"
         "The spelling has been changed for consistency with similar occurrences \n"
         "in the input files.  The object is still named group internal to \n"
         "collection.c, because it does mean something slightly different from the \n"
         "groups objects in other locations.\n\n"
         "Simply change your keyword to from 'group' to 'groups' in the COLLECTION object.\n\n"
         );

   abortAll(-1);
}




/* Local Variables: */
/* tab-width: 3 */
/* End: */


/* Local Variables: */
/* tab-width: 3 */
/* End: */
