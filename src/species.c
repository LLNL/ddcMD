#include "species.h"
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "three_algebra.h"
#include "object.h"
#include "error.h"
#include "ddcMalloc.h"
#define MAX(A,B) ((A) > (B) ? (A) : (B))
static SPECIES **specieslist = NULL;
static int namelistsize=0; 
static int typelistsize=0; 
static char **namelist=NULL; 
static char **typelist=NULL; 
static int nspecies=0; 
SPECIES *species_init(void *parent, char *name)
{
	SPECIES *sp;
	char *type;
	sp = (SPECIES *) object_initialize(name, "SPECIES", sizeof(SPECIES));
	sp->parent = parent; 
	sp->nMember=0; 
   sp->energyInfo.virialCorrection = NULL; 
	object_get((OBJECT *) sp, "type", &type, STRING, 1, "");
	if (strncmp(type, "ATOM", strlen("ATOM")) == 0)
	{
		ATOMTYPE_PARMS *parms;
		parms = ddcCalloc(1, sizeof(ATOMTYPE_PARMS));
		sp->natoms = 1;
		sp->itype = ATOMTYPE;
		sp->type = strdup(type);
 		sp->index = nspecies++;
		object_get((OBJECT *) sp, "mass", &parms->mass, WITH_UNITS, 1, "1.0","m",NULL);
		object_get((OBJECT *) sp, "charge", &parms->charge, WITH_UNITS, 1, "0.0","i*t",NULL);
		sp->parm = (void *)parms;  
	}
	return sp;
}

void resetStaticSpeciesVars()
{
   namelistsize=0;
   typelistsize=0; 
   namelist=NULL;
   typelist=NULL;
   nspecies=0;

}

void species_get(SPECIES*species, int get, void *ptr)
{
	int i, j, k, max;
	switch (get)
	{
	case NATOMS:
		*(int *)ptr = species->natoms;
		return;
	case SPECIESTYPE:
		*((char **)ptr) = species->type;
		return;
	case SPECIESITYPE:
		*((int *)ptr) = species->itype;
		return;
	case NSPECIES:
	        *(int*)ptr = nspecies;
		return;
	case SPECIESMAXNAMELENGTH:
		i = max = 0;
		while (specieslist[i] != NULL)
		{
			max = MAX(max, (int)strlen(specieslist[i]->name));
			i++;
		}
		*(int *)ptr = max;
		return;
	case SPECIESMAXTYPELENGTH:
		i = max = 0;
		while (specieslist[i] != NULL)
		{
			max = MAX(max, (int)strlen(specieslist[i]->type));
			i++;
		}
		*(int *)ptr = max;
		return;
	case NAMELIST:
		i = 0;
		while (specieslist[i] != NULL) i++;
		i++;
		if (i > namelistsize || namelist==NULL) 
		{
			namelistsize=i; 
			namelist = ddcRealloc(namelist,namelistsize*sizeof(char *));
		}
		i=0; 
		while (specieslist[i] != NULL) 
		{
			namelist[i]=specieslist[i]->name;
			i++;
		}
		namelist[i] = NULL; 
		*(char ***)ptr = namelist; 
		return;
	case TYPELIST:
		i=k=0;
		while (specieslist[k] != NULL) 
		{
			char *type; 
			type=specieslist[k++]->type;
			for (j=0;j<i;j++) if (strcmp(type,typelist[j])==0) break;
			if (j==i || i==0) 
			{
				if (i+2 > typelistsize) 
				{
					typelistsize=i+2; 
					typelist = ddcRealloc(typelist,typelistsize*sizeof(char *));
				}
				typelist[i++]=type;
			}
		}
		typelist[i] = NULL; 
		*(char ***)ptr = typelist; 
		return;
	  case SPECIESLIST:
	        *((SPECIES***)ptr) = specieslist;
	   return;
	  default:
	   assert(0);
	   break;
	}
}

int species_put(SPECIES*species, int put, void **ptr)
{
	switch (put)
	{
	case SPECIESLIST:
		specieslist = (SPECIES **) ptr;
		return 0;
	default:
		return 0;
	}
}

SPECIES *species_find(SPECIES ** species, char *name)
{
	int i;
	i = 0;
	if (species == NULL) species = specieslist;
	while (species[i] != NULL)
	{
		if (strcmp(species[i]->name, name) == 0) return species[i];
		i++;
	}
	return NULL;
}

SPECIES *species_by_index(SPECIES ** species, int index)
{
	int i;
	i = 0;
	if (species == NULL) species = specieslist;
	if (species[index]->index == index) return species[index];
	while (species[i] != NULL)
	{
		if (species[i]->index == index) return species[i];
		i++;
	}
	return NULL;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
