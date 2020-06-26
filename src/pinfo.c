#include "pinfo.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "error.h"
#include "species.h"
#include "group.h"
#include "ddcMalloc.h"

/** Returns 1 if the string is a member of the list, 0 otherwise. */
static unsigned find(char* string, char** list, unsigned nList);
static void pinfoInit(PINFO_CODEC*);


PINFO_CODEC* pinfoEncodeInit(void)
{
   PINFO_CODEC* p = ddcMalloc(sizeof(PINFO_CODEC));
   pinfoInit(p);

   GROUP** groupList;
   group_get(NULL, GROUPLIST, (void*)&groupList);
   int nGroups = 0;
   while (groupList[nGroups] != NULL)
      nGroups++;
   
   p->_gNames = ddcMalloc(nGroups * sizeof(char*));
   p->_gMap = ddcMalloc(nGroups * sizeof(unsigned));
   
   p->_nGroups=0;
   for (int ii=0; ii<nGroups; ++ii)
   {
      if (find(groupList[ii]->name, p->_gNames, p->_nGroups))
	 continue;
      p->_gNames[p->_nGroups] = strdup(groupList[ii]->name);
      assert(groupList[ii]->index < nGroups);
      p->_gMap[groupList[ii]->index] = p->_nGroups; 
      ++p->_nGroups;
   }
   
   SPECIES** speciesList;
   species_get(NULL, SPECIESLIST, (void*)&speciesList);
   int nSpecies = 0;
   while (speciesList[nSpecies] != NULL)
      ++nSpecies;

   p->_sNames = ddcMalloc(nSpecies * sizeof(char*));
   p->_sMap = ddcMalloc(nSpecies * sizeof(unsigned));
   
   p->_nSpecies = 0;
   for (int ii=0; ii<nSpecies; ++ii)
   {
      if (find(speciesList[ii]->name, p->_sNames, p->_nSpecies))
	 continue;
      p->_sNames[p->_nSpecies] = strdup(speciesList[ii]->name);
      assert(speciesList[ii]->index < nSpecies);
      p->_sMap[speciesList[ii]->index] = p->_nSpecies;
      ++p->_nSpecies;
   }
   
   p->_tNames = ddcMalloc(nSpecies * sizeof(char*));
   p->_tMap = ddcMalloc(nSpecies * sizeof(unsigned));
   
   p->_nTypes = 0;
   for (int ii=0; ii<nSpecies; ++ii)
   {
      if (find(speciesList[ii]->type, p->_tNames, p->_nTypes))
	 continue;
      p->_tNames[p->_nTypes] = strdup(speciesList[ii]->type);
      assert((int)speciesList[ii]->itype < nSpecies);
      p->_tMap[speciesList[ii]->itype] = p->_nTypes;
      ++p->_nTypes;
   }

   return p;
}

PINFO_CODEC* pinfoDecodeInit(OBJECT* obj)
{
   PINFO_CODEC* p = ddcMalloc(sizeof(PINFO_CODEC));
   pinfoInit(p);
   
   p->_nGroups = object_getv(obj, "groups", (void*)&p->_gNames, STRING,IGNORE_IF_NOT_FOUND);
   p->_nSpecies = object_getv(obj, "species", (void*)&p->_sNames, STRING,IGNORE_IF_NOT_FOUND);
   p->_nTypes = object_getv(obj, "types", (void*)&p->_tNames, STRING,IGNORE_IF_NOT_FOUND);
   
   p->_groupMap = ddcMalloc(p->_nGroups*sizeof(GROUP*));
   for (unsigned ii=0; ii<p->_nGroups; ++ii)
      p->_groupMap[ii] = group_find(NULL, p->_gNames[ii]);

   p->_speciesMap = ddcMalloc(p->_nSpecies*sizeof(SPECIES*));
   for (unsigned ii=0; ii<p->_nSpecies; ++ii)
      p->_speciesMap[ii] = species_find(NULL, p->_sNames[ii]);

return p;
}


void pinfoFree(PINFO_CODEC* p)
{
   if (!p)
      return;
   for (unsigned ii=0; ii<p->_nGroups; ++ii)
      ddcFree(p->_gNames[ii]);
   for (unsigned ii=0; ii<p->_nSpecies; ++ii)
      ddcFree(p->_sNames[ii]);
   for (unsigned ii=0; ii<p->_nTypes; ++ii)
      ddcFree(p->_tNames[ii]);
   ddcFree(p->_gNames);
   ddcFree(p->_sNames);
   ddcFree(p->_tNames);
   ddcFree(p->_gMap);
   ddcFree(p->_sMap);
   ddcFree(p->_tMap);
   ddcFree(p->_groupMap);
   ddcFree(p->_speciesMap);
   ddcFree(p);
}

LONG64 pinfoEncode(GROUP* group, SPECIES* species, PINFO_CODEC* p)
{
   unsigned iGroup   = p->_gMap[group->index];
   unsigned jSpecies = p->_sMap[species->index];
   unsigned kSType   = p->_tMap[species->itype];
   
   return  iGroup + (jSpecies+kSType)*p->_nGroups + kSType*p->_nSpecies;
}

/** Returns -1 if input (in) is outside the range that can be successfully decoded. */
int
pinfoDecode(unsigned in, GROUP** group, SPECIES** species, char** type, PINFO_CODEC* p)
{
   unsigned iMax = p->_nGroups * p->_nSpecies * p->_nTypes;
   if (in >= iMax)
      return -1;

   unsigned iGroup   = in % p->_nGroups;
   in /= p->_nGroups;
   unsigned iSpecies = in % p->_nSpecies;
   in /= p->_nSpecies;
   unsigned iType = in;

   *group   = p->_groupMap[iGroup];
   *species = p->_speciesMap[iSpecies];
   *type = p->_tNames[iType];
   return 0;
}

LONG64 pinfoMaxIndex(PINFO_CODEC* p)
{
   return p->_nGroups * p->_nSpecies * p->_nTypes;
}


/** Returns 1 if the string is a member of the list, 0 otherwise. */
unsigned find(char* string, char** list, unsigned nList)
{
   for (unsigned ii=0; ii<nList; ++ii)
      if (strcmp(string, list[ii]) == 0)
	 return 1;
   return 0;
}

void
pinfoInit(PINFO_CODEC* p)
{
   p->_nGroups=0;
   p->_nSpecies=0;
   p->_nTypes=0;
   p->_gNames=NULL;
   p->_sNames=NULL;
   p->_tNames=NULL;
   p->_gMap=NULL;
   p->_sMap=NULL;
   p->_tMap=NULL;
   p->_groupMap=NULL;
   p->_speciesMap=NULL;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
