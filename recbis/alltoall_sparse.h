/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef ALLTOALL_SPARSE__
#define ALLTOALL_SPARSE__

#include <mpi.h>
#include "integer.h"

integer alltoall_sparse(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
			integer nloc,integer nlocmax,
			integer sz,void *data,integer (*get_dest)(const void *element));

/*
Sparse version of MPI_Alltoallv(). Given on each processor a set of objects with defined
destinations, this routine ensures that all objects are delivered to their destinations.
Note that no pre-negotiation of who is sending how much data to whom is required. Each
processor will perform O(lop(P)**q) sends and receives, where P is the total number of
participating processors (P = (pid1-pid0)/pstride). If the equalize() routine is run all
the time (see source) then q = 2, if it run only infrequently, q = 1.

The problem is solved by a group of processors, with MPI ranks:

   [pid0 , pid0+pstride , pid0+2*pstride , . . . , pid1).

Arguments:
   pid0    - - lowest processor MPI rank in group
   pid1    - - 1+higest processor MPI rank in group
   pid     - - MPI rank of this processor   
   pstride - - rank spacing between processor in group.
   comm    - - MPI communicator used for communication and which
               ranks pertain to

   nloc    - - number of elements on this processor on entry
   nlocmax - - number of elements the input array can hold, must be large enough
   sz      - - size of each element (e.g. from sizeof())
   data    - - pointer to elements
   get_dest- - given a pointer to an element, return its destination rank

Return value:
   Returns the current number of elements in `data', i.e. the number of elements that
   had this processor as destination.
*/


#endif
