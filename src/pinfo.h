#ifndef PINFO_H
#define PINFO_H

#include "gid.h"
#include "object.h"
#include "group.h"
#include "species.h"

typedef struct pinfoCodec_st
{
   unsigned  _nGroups;
   unsigned  _nSpecies;
   unsigned  _nTypes;
   char**    _gNames;
   char**    _sNames;
   char**    _tNames;
   unsigned* _gMap;
   unsigned* _sMap;
   unsigned* _tMap;
   GROUP**   _groupMap;
   SPECIES** _speciesMap;
} PINFO_CODEC;


PINFO_CODEC* pinfoEncodeInit(void);
PINFO_CODEC* pinfoDecodeInit(OBJECT*);
void pinfoFree(PINFO_CODEC*);
LONG64 pinfoEncode(GROUP*, SPECIES*, PINFO_CODEC*);
int pinfoDecode(unsigned in, GROUP**, SPECIES**, char** type, PINFO_CODEC*);
LONG64 pinfoMaxIndex(PINFO_CODEC*);

#endif // #ifndef PINFO_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
