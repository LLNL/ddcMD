/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef EQUALIZE__
#define EQUALIZE__

#include <mpi.h>

#include "integer.h"

integer equalize(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc_p[1],integer nlocmax,integer sz,void *data_p);

/*
Equalization of data elements across processors. After this call, the maximum and minimum
number of elements among the processors will differ by at most one.


The problem is solved by a group of processors, with MPI ranks:

   [pid0 , pid0+pstride , pid0+2*pstride , . . . , pid1).

Arguments:
   pid0    -- lowest processor MPI rank in group
   pid1    -- 1+higest processor MPI rank in group
   pid     -- MPI rank of this processor   
   pstride -- rank spacing between processor in group.
   comm    -- MPI communicator used for communication and which
              ranks pertain to
   
   nloc_p  -- Pointer to integer. That integer contains on entry the number of elements
               on this processor. On return, it contains the new number of elements on
               this processor.
   nlocmax -- The number of elements the input array can hold (amount of allocated
              elements). Must be sufficiently large for subroutine to succeed.
   sz      -- The size of each element (e.g. from sizeof()).
   data    -- Pointer to elements

Return value:
   The total number of elements (aggregated over all processors) is returned.

The algorithm is recursive with about log_2 ((pid1-pid0)/pstride) steps. If the
data is relatively evenly distributed to begin with, very little data is exchanged.
*/

#endif
