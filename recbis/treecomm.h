/*
   This file is part of RECBIS:
   A library for parallel computation of medians and for
   parallel data distribution using recursive bisection.

   Written by Tomas Oppelstrup at the Lawrence Livermore
   National Laboratory.
*/

#ifndef TREECOMM__
#define TREECOMM__

#include <mpi.h>

#include "integer.h"

typedef void (*accum_fun)(void *x,const void *dx);


/*
These are a set of communication routines that operate on a strided subset of ranks in
an MPI communicator. The communication is among the group of processors with MPI ranks:

   [pid0 , pid0+pstride , pid0+2*pstride , . . . , pid1).

Common argument names and meanings:
   pid0    - - lowest processor MPI rank in group
   pid1    - - 1+higest processor MPI rank in group
   pid     - - MPI rank of this processor   
   pstride - - rank spacing between processor in group.
   comm    - - MPI communicator used for communication and which
               ranks pertain to

   nloc,ndata - - number of elements on this processor on entry
   sz         - - size of each element (e.g. from sizeof())
   data       - - pointer to elements

Description of individual routines

tree_gather_inc:
   Gather data on first processor in group (pid0). Each processor can contribute
   up to nloc_max elements (must be same on all participating processors. nloc
   is the number of elements provided by this processor. On exit, processor with
   rank pid0 will have in data all elements from all processors. data needs to be
   able to hold at least (pid1-pid0)/pstride*nloc_max elements. Returns the number
   of elements in data after gathering. On pid0 this is the total number of
   elements provided by all processors. On other processors this will be a
   possibly smaller number.

tree_gather:
   Gather data on first processor in group (pid0). Each processor contributes
   nloc elements (must be same on all participating processors). On return,
   processor with rank pid0 will have in data all elements provided by all
   processors.

tree_broadcast:
   Elements on pid0 are broadcast to all participating processors. When the call
   returns, that processor has the same content in data as initially on pid0.

tree_reduce:
   Like MPI_Reduce. acc is a function that given pointers to two elements, accumulates
   the second element into the first, e.g. for summation of double precision numbers:

       void acc_double(void *a,const void *b) {
         ((double *) a)[0] += ((const double *) b)[0];
       }
      
   When pid0's call returns, element j of its data array contains the accumulation of
   element j on all processors. For other processors, their data arrays contains partial
   accumulations (e.g. like an inclusive scan).
   
tree_allreduce:
   Like MPI_Allreduce. Same as tree_reduce() above, except all processors contain
   the complete reduction on return.

tree_barrier:
   Like MPI_Barrier. Blocks until all participating processor have entered the barrier.
*/

integer tree_gather_inc(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
			integer nloc_max,integer nloc,integer sz,void *data);

void tree_gather(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer nloc,integer sz,void *data);

void tree_broadcast(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer ndata,integer sz,void *data);

void tree_reduce(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		 integer ndata,integer sz,void *data,accum_fun acc);

void tree_allreduce(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm,
		    integer ndata,integer sz,void *data,accum_fun acc);

void tree_barrier(integer pid0,integer pid1,integer pid,integer pstride,MPI_Comm comm);
#endif
