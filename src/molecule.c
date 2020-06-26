#include "molecule.h"
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "ddcMalloc.h"
#include "system.h"
#include "bioGid.h"
#include "heap.h"
#include "expandbuffer.h"
#include "ddc.h"
#include "preduce.h"
#include "mpiUtils.h"
static int nMolecule=0; 
MOLECULETYPE *moleculeInit(void *parent, char *name);
MOLECULECLASS *moleculeClassInit(void *parent, char *name)
{
	MOLECULECLASS *moleculeClass = (MOLECULECLASS *) object_initialize(name, "MOLECULECLASS", sizeof(MOLECULECLASS));
	moleculeClass->parent = parent; 
   moleculeClass->nMoleculeTypes =0; 
   moleculeClass->moleculeTypes = NULL; 
   moleculeClass->speciesToMoleculeOffset = NULL; 
   moleculeClass->nSpecies = 0 ;
   moleculeClass->species = NULL; 
   moleculeClass->list = NULL; 
   moleculeClass->molecule = NULL; 
   char *type;  
   object_get((OBJECT *) moleculeClass, "type", &type, STRING, 1, "NORMAL");
   moleculeClass->name = ddcCalloc(strlen(name) + 1, sizeof(char));
   moleculeClass->type = ddcCalloc(strlen(type) + 1, sizeof(char));
   strcpy(moleculeClass->name, name);
   strcpy(moleculeClass->type, type);

   char **moleculeNames=NULL; 
   moleculeClass->nMoleculeTypes=object_getv((OBJECT *) moleculeClass, "molecules",   (void *)&moleculeNames,   STRING, ABORT_IF_NOT_FOUND);

   moleculeClass->moleculeTypes = (MOLECULETYPE **)ddcMalloc(moleculeClass->nMoleculeTypes*sizeof(MOLECULETYPE*)); 
   int nSpecies =0; 
   for (int i=0;i<moleculeClass->nMoleculeTypes;i++) 
   { 
      OBJECT *obj = object_find(moleculeNames[i],"MOLECULE");
      int size = object_keywordSize(obj,"species"); 
      nSpecies+= size; 
   }
   moleculeClass->nSpecies=nSpecies; 
   moleculeClass->species=(SPECIES **)ddcMalloc((nSpecies+1)*sizeof(SPECIES *)); 
   moleculeClass->speciesBuffer = moleculeClass->species;
   moleculeClass->species[moleculeClass->nSpecies] = NULL;  //end of list marker. 
   moleculeClass->speciesIndexToMoleculeIndex= (int *)ddcMalloc((nSpecies)*sizeof(int)); 


   for (int i=0;i<moleculeClass->nMoleculeTypes;i++) 
   {
      MOLECULETYPE *molecule =moleculeClass->moleculeTypes[i]=moleculeInit(moleculeClass,moleculeNames[i]);
      for (int j=0;j<molecule->nSpecies;j++) moleculeClass->speciesIndexToMoleculeIndex[molecule->species[j]->index] = molecule->index;
      moleculeClass->speciesBuffer += molecule->nSpecies; //Advance the buffer for next molecule; 
   }
   ddcFree(moleculeNames); 


	return moleculeClass;
}
void moleculeDelta(MOLECULECLASS *moleculeClass, STATE *state, THREE_VECTOR *delta)
{
   MOLECULE *molecule =  moleculeClass->molecule; 
   int nMolecules =  moleculeClass->nMolecules; 
   for (int k=0;k<nMolecules;k++)
   {
      MOLECULETYPE  *type = molecule[k].type; 
      int a = type->ownershipSpeciesOffset; 
      int i0 = molecule[k].list[a];
      THREE_VECTOR r0 = {state->rx[i0],state->ry[i0],state->rz[i0]};
      double M=0.0; 
      THREE_VECTOR R=vzero; 
      THREE_VECTOR d[type->nSpecies]; 
      for (int a=0;a<type->nSpecies;a++)
      {
         int i = molecule[k].list[a];
         double mass = ((ATOMTYPE_PARMS *)(state->species[i]->parm))->mass;
         VSET(d[a],state->rx[i]-r0.x,state->ry[i]-r0.y,state->rz[i]-r0.z); 
         nearestImage(&d[a].x,&d[a].y,&d[a].z); 
         VSVOP(R,+=,mass,*,d[a]);
         M += mass; 
      }
      VSCALE(R,1.0/M); 
      for (int a=0;a<type->nSpecies;a++)
      {
         VOP1(d[a],-=,R);
         int i = molecule[k].list[a];
         delta[i] = d[a]; 
      }
   }
}
void testSpeciestoMoleculeOffset(SYSTEM *sys)
{
   STATE *state = sys->collection->state; 
   gid_type *label=state->label; 
   MOLECULECLASS *moleculeClass = sys->moleculeClass;
   for (unsigned i=0;i<sys->nlocal;i++)
   {
      int iMoleIndex=moleculeClass->speciesToMoleculeOffset[i];
      if (iMoleIndex == -1) continue; 
      //unsigned atmI =(unsigned)(label[i]&atmgrpMask);
      //gid_type moleI= label[i]& molMask; 
      //int speciesIndex = state->species[i]->index;
      MOLECULE molecule = moleculeClass->molecule[iMoleIndex]; 
      long long unsigned gidM =molecule.gid;
      //char *moleName=molecule.type->name;
      //int typeIndex=molecule.type->index;
      /*if (moleI>>32 != gidM)*/ fprintf(ddcfile,"%d: %d %lu %llx %s\n",getRank(0),i,label[i]>>32,gidM,state->species[i]->name); 
   }
}
void moleculeScanState(MOLECULECLASS *moleculeClass, STATE *state)
{
   struct { unsigned gid; int  typeIndex;} *map; 
   unsigned mapBlk; 
   map = heapGet(&mapBlk); 
   int mapCapacity = heapBlockSize(mapBlk)/sizeof(*map); 

   int nMoleculeTypes = moleculeClass->nMoleculeTypes;
   MOLECULETYPE **moleculeTypes= moleculeClass->moleculeTypes;

   for (int i=0;i<nMoleculeTypes;i++) moleculeTypes[i]->nMembers=0; 

// count and make a list of number of multi-species molecules. 
   unsigned gidMax=0;
   int nn =0; 
   int nMoleculeGlobal=0; 
   for (int i=0;i<state->nlocal;i++)
   {
      assert(nn+1 <= mapCapacity); 
      unsigned gid= (state->label[i]>>32); 
      SPECIES *species = state->species[i]; 

      int moleculeTypeIndex = moleculeClass->speciesIndexToMoleculeIndex[species->index];
      MOLECULETYPE *type = moleculeTypes[moleculeTypeIndex]; 
      if (species == type->ownershipSpecies)
      {
          nMoleculeGlobal++;
         if(type->nSpecies > 1 ) 
         {
            moleculeTypes[moleculeTypeIndex]->nMembers++;  
            map[nn].gid = gid; 
            map[nn].typeIndex = moleculeTypeIndex; 
            if (gid > gidMax) gidMax=gid; 
            nn++; 
         }
      }
   }
   MPI_Allreduce(&nMoleculeGlobal,&moleculeClass->nMoleculesGlobal,1,MPI_INT,MPI_SUM,COMM_LOCAL);

   heapEndBlock(mapBlk, nn*sizeof(*map));


   int nMolecules=0; 
   for (int i=0;i<moleculeClass->nMoleculeTypes;i++) nMolecules +=moleculeTypes[i]->nMembers; 
   assert(nn == nMolecules); 
   moleculeClass->nMolecules = nMolecules; 
   unsigned moleculeGidToIndexBlk; 
   int *moleculeGidToIndex=heapGet(&moleculeGidToIndexBlk); 
   assert( gidMax < heapBlockSize(moleculeGidToIndexBlk)/sizeof(*moleculeGidToIndex)); 
   heapEndBlock(moleculeGidToIndexBlk, (gidMax+1)*sizeof(*moleculeGidToIndex));

   int offset[nMoleculeTypes+1];   
   int listOffset[nMoleculeTypes+1]; 
   offset[0]=0; 
   listOffset[0]=0; 
   for (int i=1;i<=nMoleculeTypes;i++) 
   {
      int cnt = moleculeTypes[i-1]->nMembers;
      offset[i] = offset[i-1] + cnt;  
      listOffset[i] = listOffset[i-1] + cnt*moleculeTypes[i-1]->nSpecies;
   }
   int listSize = listOffset[nMoleculeTypes];

   // list hold the species offset in state
   int *list =moleculeClass->list=ExpandBuffers((void*) moleculeClass->list, sizeof(int), listSize, 500, LOCATION("moleculeScanState"), "moleculeClass->list");
   // speciesToMoleculeOffset is a map from species offset in state to molecule offset in molecule; 
   moleculeClass->speciesToMoleculeOffset=ExpandBuffers((void*) moleculeClass->speciesToMoleculeOffset, sizeof(int), state->nion, 500, LOCATION("moleculeScanState"), "moleculeClass->speciesToMoleculeOffset");
   //  molecule is a molecularType ordered list of  multispecies moleculues. 
   moleculeClass->molecule=ExpandBuffers((void*) moleculeClass->molecule, sizeof(MOLECULE), nMolecules, 500, LOCATION("moleculeScanState"), "moleculeClass->molecules");
   MOLECULE *molecules=moleculeClass->molecule;

   for (int i =0;i<listSize;i++) list[i]=-1; 
   for (unsigned i=0;i<=gidMax;i++)  moleculeGidToIndex[i]=-1;  

   // Create an order list of molecules
   for (int i = 0; i<  nn; i++) 
   {
      unsigned gid = map[i].gid;
      int typeIndex = map[i].typeIndex;
      int index = offset[typeIndex]; 
      moleculeGidToIndex[gid] = index; 
      molecules[index].gid = gid;
      molecules[index].type = moleculeTypes[typeIndex];
      molecules[index].list = list+listOffset[typeIndex]; 

      offset[typeIndex]++; 
      listOffset[typeIndex]+=moleculeTypes[typeIndex]->nSpecies; 
   }

   // populate the atom list of each molecule. 
   for (int i=0;i<state->nion;i++)
   {
      unsigned gid= (state->label[i]>>32); 
      int index = -1;
      if (gid<=gidMax) index = moleculeGidToIndex[gid]; 
      moleculeClass->speciesToMoleculeOffset[i]=index; 
      if (index == -1 ) continue ; 
      unsigned atom= (state->label[i])&atmgrpMask; 
      molecules[index].list[atom]=i;
      moleculeClass->speciesToMoleculeOffset[i]=index; 
   } 
   // Clean up temporary storage on the heap
   heapFree(moleculeGidToIndexBlk); 
   heapFree(mapBlk); 
}
MOLECULETYPE *moleculeInit(void *parent, char *name)
{
   MOLECULECLASS *moleculeClass = (MOLECULECLASS *)parent; 
   MOLECULETYPE *molecule = (MOLECULETYPE *) object_initialize(name, "MOLECULE", sizeof(MOLECULETYPE));
   char **speciesNames=NULL; 
   molecule->nSpecies=object_getv((OBJECT *) molecule, "species",   (void *)&speciesNames,   STRING, ABORT_IF_NOT_FOUND);
   molecule->species = moleculeClass->speciesBuffer;
   char *ownershipSpecies;
   object_get((OBJECT *)molecule ,"ownershipSpecies", &ownershipSpecies, STRING, 1, "$NONE$");
   molecule->ownershipSpecies=NULL; 
   for (int i=0;i<molecule->nSpecies;i++)
   {
      molecule->species[i] = species_init(system,speciesNames[i]);
      moleculeClass->speciesIndexToMoleculeIndex[molecule->species[i]->index] = molecule->index;
      if (strcmp(speciesNames[i],ownershipSpecies)==0) 
      {
         molecule->ownershipSpecies=molecule->species[i]; 
         molecule->ownershipSpeciesOffset=i;
      }
   }
   if (molecule->ownershipSpecies==NULL)   
   {
      molecule->ownershipSpecies=molecule->species[0];
      molecule->ownershipSpeciesOffset=0;
   }

   molecule->nMembers=0; 
   molecule->index = nMolecule++;
   ddcFree(speciesNames); 
   return molecule; 
}



/* Local Variables: */
/* tab-width: 3 */
/* End: */
