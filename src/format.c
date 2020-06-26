#include <stdio.h>
#include <inttypes.h>
#include <strings.h>
static char _loopFormat[16];
static char _gidFormat[16];
static int _loopFormatSize;
void loopFormatInit(int nDigits) 
{ 
   _loopFormatSize=nDigits;
   sprintf(_loopFormat,"%%%d.%d"PRIu64,nDigits,nDigits); 
}
void gidFormatInit(char *base) 
{ 
   if (strcasecmp(base,"decimal")==0) sprintf(_gidFormat,"%s","%12.12"PRIu64); 
   if (strcasecmp(base,"hex")==0)     sprintf(_gidFormat,"%s","%16.16"PRIx64); 
}
char *loopFormat() { return _loopFormat; }
int loopFormatSize() { return _loopFormatSize;}

char *gidFormat()  { return _gidFormat; }
