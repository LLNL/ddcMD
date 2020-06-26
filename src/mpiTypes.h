#ifndef MPI_TYPES_H
#define MPI_TYPES_H

#include <mpi.h>

/** Self initializing functions to construct various MPI_Datatypes.
 *  These functions automatically construct the corresponding type the
 *  first time they are called.  This eliminates the need for explicit
 *  initialization.
 * 
 *  You might be tempted to declare a function like
 *  threeVector_MPIType() in three_algebra.h instead of here.  However,
 *  doing so would require you to include mpi.h in three_algebra.h and
 *  you really don't want to couple everything that needs a THREE_VECTOR
 *  to mpi.  On the other hand, declaring all of these functions here
 *  does not require any additional header files (types such as
 *  PARTICLE are needed in mpiTypes.c only). */

MPI_Datatype domainx_MPIType(void);
MPI_Datatype threeVector_MPIType(void);
MPI_Datatype threeInt_MPIType(void);
MPI_Datatype particle_MPIType(void);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
