#ifndef MPI_ALGORITHM_H
#define MPI_ALGORITHM_H

#include <mpi.h>
#include "gid.h" // LONG64

/** Compute the maximum value of the distributed array data.  nData is
 *  the number of elements of data owned by a given task (this need not
 *  be the same on all tasks). */
LONG64 mpiMaxVal(LONG64* data, unsigned nData,MPI_Comm);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
