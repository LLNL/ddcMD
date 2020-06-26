/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef BALANCE__
#define BALANCE__

#include <mpi.h>
#include "integer.h"
#include "genptr.h"

integer balance2(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc,integer nlocmax,integer sz,void *data,cmp_fun cmp,
		 void *pivot);

integer balance(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		integer nloc,integer nlocmax,integer sz,void *data,cmp_fun cmp);

/*
Compute a splitting point so half the data elements are to to the left and half to the right.
Data is exchanged so that the first half of the processors hold the left half, and the second
half of the processors hold the right half. If the number of processors is an odd number, the
splitting point is adjusted so that each of the two sets of processors has a number of elements
proportional to the number of processors in that set. The only difference between balance() and
balance2() is that balance2() returns the chosen splitting point.

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
   cmp     - - comparison function, so that cmp(a,b) returns <0, 0, >0 depending on
               whether a<b, a=b, or a>b, respectively.
   pivot   - - pointer to element (not in data!). On return, that element will contain
               (a copy of) the pivot element (splitting point).

Return value:
   Returns the number of elements currently on this processor (in the data array).
*/

#endif
