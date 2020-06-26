#ifndef NEXTFILE_H
#define NEXTFILE_H
#include "gid.h"
#include "pio.h"
typedef struct nextfile_parms_st
{
	SIGNED64 loopindex;
	int nfiles;
	char **files;
   int *filetype;
   PFILE *pfile;
	char *currentfile;
} NEXTFILE_PARMS;
#endif
