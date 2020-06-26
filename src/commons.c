#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/stat.h>
#include "ddc.h"

extern FILE *ddcfile;
FILE *getddcfile()
{
	return ddcfile;
}

void commons_init(void)
{
	int domainId;
	char filename[256];
	struct stat statbuf;
	MPI_Comm_rank(MPI_COMM_WORLD, &domainId);
	ddcfile = NULL; 
   int rc=-1;
   if (domainId == 0) rc = stat("ddcfiles", &statbuf);
   MPI_Bcast(&rc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rc == 0)
	{
		sprintf(filename, "ddcfiles/ddcfile%6.6d", domainId);
		ddcfile = fopen(filename, "w");  
	}
}

/* Local Variables: */
/* tab-width: 3 */
/* End: */
