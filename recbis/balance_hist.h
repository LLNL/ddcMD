/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef BALANCE_HIST__
#define BALANCE_HIST__

#include "treecomm.h"
#include "bisection_data.h"

integer balance_hist(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		     integer nloc,integer nlocmax,bisection_data data[],
		     integer idir,double x0,double x1,double volume_cost,
		     double xcut_out[1]);
/*
Compute a splitting point so half the data elements are to to the left and half to the right.
Data is exchanged so that the first half of the processors hold the left half, and the second
half of the processors hold the right half. If the number of processors is an odd number, the
splitting point is adjusted so that each of the two sets of processors has a number of elements
proportional to the number of processors in that set. This is akin to balance() and balance2(),
with the difference that data are explicit 3D points, and each point as a weight. The splitting
is computed along one of the X-, Y-, and Z- Cartesian directions.

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
   idir    - - Cartsian direction along which to split (0 = x, 1 = Y, 2 = Z).
   x0      - - leftmost value (interval beginning) in `idir'-direction
   x1      - - rightmost value (interval end) in `idir'-direction
   volume_cost - - additional cost per unit length along `idir'-direction
   xcut_out    - - pointer to floating point number to which splitting point
                   will be stored on return.

Return value:
   Returns the number of elements currently on this processor (in the data array).
*/

#endif
