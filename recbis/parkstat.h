/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef PARKSTAT__
#define  PARKSTAT__

#include <mpi.h>
#include "integer.h"

/* This defines the function pointer type cmp_fun */
#include "genptr.h"

integer parkstat(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer k,integer nloc,integer sz,void *vloc,void *pivot,cmp_fun cmp);

/*
Parallel computation of the k-statistic of a data set, meaning if the
data set were to be sorted in increasing order, return the k'th element.
See kstat() in kstat.h for a serial implementation.

The problem is solved by a group of processors, with MPI ranks:

   [pid0 , pid0+pstride , pid0+2*pstride , ... , pid1).

Arguments:
   pid0    -- lowest processor MPI rank in group
   pid1    -- 1+higest processor MPI rank in group
   pid     -- MPI rank of this processor   
   pstride -- rank spacing between processor in group.
   comm    -- MPI communicator used for communication and which
              ranks pertain to

   k     -- which element to return
   nloc  -- number of elements in input array on this processor
   sz    -- size of each element (e.g. as from sizeof() operator)
   vloc  -- pointer to elements on this processor
   cmp   -- function pointer to comparison function
   pivot -- pointer to element, which on return will be set to k-statsitic.

Return value:
   The number of pivot elements is the returned, i.e. the total number of
   elements across all processors that compare equal to the returned
   pivot element (returned in pivot).

This routine is estimated to run in time proportional to n/p + log(n)*log(p),
where nlmax is the maximum number of elements on any one processor, and p is
the number of processors. It is assumed that cmp takes constant time.
*/

#endif
