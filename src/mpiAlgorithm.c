#include "mpiAlgorithm.h"
#include <mpi.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))

LONG64 mpiMaxVal(LONG64* data, unsigned nData, MPI_Comm comm)
{
   LONG64 localMax=0;
   LONG64 globalMax=0;
   for (unsigned ii=0; ii<nData; ++ii)
      localMax = MAX(localMax, data[ii]);
   
   MPI_Allreduce(&localMax, &globalMax, 1, MPI_LONG_LONG, MPI_MAX, comm);
   
   return globalMax;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
