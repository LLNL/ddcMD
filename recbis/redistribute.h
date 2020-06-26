/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef REDISTRIBUTE__
#define REDISTRIBUTE__

#include <mpi.h>

#include "integer.h"
#include "genptr.h"

integer redistribute(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		     integer nloc,integer nlocmax,integer sz,void *data,
		     int ndir,int idir,cmp_fun cmp_coord[/*ndir*/]);

/*
Perform a recursive bisection partitioning of K dimensional point data, so that on return,
each processor has as close as possible the same number of points confined in a rectangular
brick.  The splitting direction cycles over the K dimensions.

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
   ndir    - - dimensionality of points (K in description above)
   idir    - - first dimension to split along
   cmp_coord  - - Array of `ndir' comparison functions. The cmp_coord[d] should
                  compare dimension d of the input points.

Return value:
   Returns the number of points on this processor (in data array) after the computation.
*/

#endif
